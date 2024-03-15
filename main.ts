/**
 * An extension providing a compass-heading for a buggy located anywhere on the globe 
 * (except at the magnetic poles!), with any mounting orientation for the microbit.
 */

//% color=#6080e0 weight=40 icon="\uf14e" block="Heading" 
namespace heading {
    enum Dim { // ...for brevity
        X = Dimension.X,
        Y = Dimension.Y,
        Z = Dimension.Z
    }

    // GLOBALS
    const MarginalField = 30 // minimum acceptable field-strength for magnetometer readings
    const MinPeakSeparation = 500 // sanity check to ignore close peaks due to wobbling signal
    // (still permits maximum spin-rate of 120 RPM, or 2 rotations a second!) 
    const Inertia = 0.95 // ratio of old to new radius readings (for inertial smoothing)
    const Window = 200 // minimum ms separation of radius extrema (~ quarter rotation time)

    let scanTimes: number[] = [] // sequence of time-stamps for scanned readings 
    let scanData: number[][]= [] // scanned sequence of [X,Y,Z] magnetometer readings  
    let uDim = -1 // the "horizontal" axis (pointing East) for transformed readings (called U)
    let vDim = -1 // the "vertical" axis (pointing North) for transformed readings (called V)
    //let wDim = -1 // the third, unused, dimension (used to assess rotation-speed)
    let uOff = 0 // the offset needed to re-centre readings along the U-axis
    let vOff = 0 // the offset needed to re-centre readings along the V-axis
    let theta = 0 // the angle to rotate readings so the projected ellipse aligns with the U & V axes
    let scale = 0 // the scaling to stretch V-axis readings from the ellipse onto a circle
    let toNorth = 0 // the angular bias to be added (so that North = 0)
    let strength = 0 // the average magnetic field-strength observed by the magnetometer
    let period = 0 // average rotation time derived from scanData[]
    // spacing of samples from scanData to test with (~ #samples covering an octant)
    let turn45 = 0 
    //let uFlip = 1 // set to -1 if uDim polarity is inverted
    //let vFlip = 1 // set to -1 if vDim polarity is inverted
    let plane: string = "**" // the projection plane we are using: "XY","YZ" or "ZX"
    let logging = true  // logging mode flag
    let testing = false  // test mode flag
    let test = 0 // selector for test sample

    // EXPORTED USER INTERFACES   

    /** 
     * Assuming the buggy is currently spinning clockwise on the spot, capture a 
     * time-stamped sequence of magnetometer readings from which to set up the compass.
     * NOTE that once scanning is complete, the heading.setNorth() function must then 
     * be called (to process this scanned data) before heading.degrees() will work.
     * @param ms scanning-time in millisecs (long enough for more than one full rotation) 
    */
    //% block="scan for (ms) $ms" 
    //% inlineInputMode=inline 
    //% ms.shadow="timePicker" 
    //% ms.defl=0 
    //% weight=90 

    export function scan(ms: number) {
        // Every ~30 ms over the specified duration (generally a couple of seconds),
        // magnetometer readings are sampled and a new [X,Y,Z] triple added to the scanData[] array.
        // A timestamp for each sample is also recorded in the scanTimes[] array.
        if (testing) {
            simulateScan()
            basic.pause(ms)
        } else {

            // (To smooth out jitter, each reading is always a rolling sum of SEVEN consecutive readings!)
            let now = input.runningTime()
            let finish = now + ms
            let index = 0
            let xRoll: number[] = []
            let yRoll: number[] = []
            let zRoll: number[] = []
            // take the first six readings...
            for (let k = 0; k < 6; k++) {
                xRoll.push(input.magneticForce(0))
                yRoll.push(input.magneticForce(1))
                zRoll.push(input.magneticForce(2))
                basic.pause(25)
            }

            let x = 0
            let y = 0
            let z = 0
            // continue cranking out rolling sums, adding a new reading and dropping the oldest
            // to each of the
            while (now < finish) {
                now = input.runningTime()
                scanTimes.push(now - 100) // the time of the middle readings (roughly)
                basic.pause(25)
                xRoll.push(input.magneticForce(0))
                x = 0
                xRoll.forEach(a => x += a) // collect x sum
                xRoll.shift()
                yRoll.push(input.magneticForce(1))
                y = 0
                yRoll.forEach(a => y += a) // collect y sum
                yRoll.shift()
                zRoll.push(input.magneticForce(2))
                z = 0
                zRoll.forEach(a => z += a) // collect z sum
                zRoll.shift()
                scanData.push([x,y,z]) // add new reading
            }
        }
        
        if (logging) {
            datalogger.deleteLog()
            datalogger.includeTimestamp(FlashLogTimeStampFormat.None)
            datalogger.setColumnTitles("index","t","x","y","z")
            for (let i = 0; i < scanTimes.length; i++) {
                datalogger.log(
                    datalogger.createCV("index", i),
                    datalogger.createCV("t", scanTimes[i]),
                    datalogger.createCV("x", scanData[i][Dim.X]),
                    datalogger.createCV("y", scanData[i][Dim.Y]),
                    datalogger.createCV("z", scanData[i][Dim.Z]))
            }
        }
    }


    /**
     * Analyse scanned data to prepare for reading compass-headings.
     * Returns either the spin RPM, or a negative error code:
     * 
     *      -1 : NOT ENOUGH SCAN DATA
     * 
     *      -2 : FIELD STRENGTH TOO WEAK
     * 
     *      -3 : NOT ENOUGH SCAN ROTATION
     */
    //% block="set North" 
    //% inlineInputMode=inline 
    //% weight=80 
    export function setNorth(): number {
    /* Analyse the scan-data to decide how best to use the magnetometer readings.

    As the buggy spins, the magnetic field-vector sweeps out a Spin-Circle on the surface of a sphere.
    In the fully general case, this projects onto the plane of each pair of orthogonal axes (XY,YZ,ZX) 
    as an Ellipse with a certain eccentricity (and maybe a tilt). We will get the best heading 
    discrimination from the "roundest", most "square-on" view (i.e. the least eccentric Ellipse).
    
    Taking readings from just those two axes, we'll need to transform points lying around the Ellipse 
    (by first rotating about the origin, and then scaling vertically to "un-squash" the Ellipse), so that 
    they lie back on the unprojected Spin-Circle, giving a relative angle that can (eventually) be offset 
    by a fixed bias to return the true heading with respect to North.
    */

        // we need at least ~3 second's worth of scanned readings...
        let nSamples = scanTimes.length
        if (nSamples < 100) {
            return -1 // "NOT ENOUGH SCAN DATA"
        }
        // Each dimension tracks a sinusoidal wave of values (generally not centred on zero).
        // The first pass finds the ranges for each axis 
        let xlo = 999
        let ylo = 999
        let zlo = 999
        let xhi = -999
        let yhi = -999
        let zhi = -999
        // To assess the range and offset, find the raw extremes.
        for (let i = 0; i < nSamples; i++) {
            xhi = Math.max(xhi, scanData[i][Dim.X])
            yhi = Math.max(yhi, scanData[i][Dim.Y])
            zhi = Math.max(zhi, scanData[i][Dim.Z])
            xlo = Math.min(xlo, scanData[i][Dim.X])
            ylo = Math.min(ylo, scanData[i][Dim.Y])
            zlo = Math.min(zlo, scanData[i][Dim.X])
        }

        // Use the mean of these extremes as normalisation offsets
        let xOff = (xhi + xlo) / 2
        let yOff = (yhi + ylo) / 2
        let zOff = (zhi + zlo) / 2

        // Now assess eccentricity by looking for the shortest and longest radii for each Ellipse.
        let xyhi = 0
        let yzhi = 0
        let zxhi = 0
        // (An initial upper limit to radius is quarter of the perimeter of its bounding-box)
        let xHalf = (xhi - xlo) / 2
        let yHalf = (yhi - ylo) / 2
        let zHalf = (zhi - zlo) / 2
        let xylo = xHalf + yHalf
        let yzlo = yHalf + zHalf
        let zxlo = zHalf + xHalf

        // square these upper bounds, which will subsequently shrink to the minimum
        xylo *= xylo
        yzlo *= yzlo
        zxlo *= zxlo
        
        // Record the angle each major radius makes (anti-clockwise from the horizontal)
        let xya = 0
        let yza = 0
        let zxa = 0

        // initialise global field-strength & the smoothed radii-squared
        strength = 0
        let xyRsq = 0
        let yzRsq = 0
        let zxRsq = 0

        // arrays for recording peak-radius sample-indices, added as we discover them
        let xyPeaks: number[] = []
        let yzPeaks: number[] = []
        let zxPeaks: number[] = []

        // initialase timestamps of most recent peaks detected
        let xyLast = scanTimes[0]
        let yzLast = scanTimes[0]
        let zxLast = scanTimes[0]

        if (logging) {
            // prepare for analysis
            datalogger.setColumnTitles("index", "x", "y", "z", "xyhi", "yzhi", "zxhi"," xya", "yza","zxa")
        }

        for (let j = 0; j < nSamples; j++) {
            // extract and normalise the next sample readings
            let stamp = scanTimes[j]
            let x = scanData[j][Dim.X] - xOff
            let y = scanData[j][Dim.Y] - yOff
            let z = scanData[j][Dim.Z] - zOff
            // By Pythagoras: radius-squared = sum of squares of coordinates
            let xsq = x * x
            let ysq = y * y
            let zsq = z * z

            // accumulate square of field-strength (a global)
            strength += xsq + ysq + zsq 

            // projection in XY plane...
            let rsq = xsq + ysq
            // in tracking the radius, we use inertial smoothing to reduce 
            // multiple detections due to minor fluctuations in readings
            xyRsq = rsq - Inertia * (rsq - xyRsq) // === xyRsq*Inertia + rsq*(1-Inertia)
            xylo = Math.min(xylo, xyRsq) // shortest so far...
            if (xyRsq > xyhi) {
                xyhi = xyRsq // longest so far...
                xya = Math.atan2(y, x) // angle (anticlockwise from X axis)
                // need to clock new peak?
                if ((stamp - xyLast) > MinPeakSeparation)
                {
                    xyPeaks.push(j)
                    xyLast = stamp
                    if (logging) {
                        datalogger.log(
                            datalogger.createCV("index", j),
                            datalogger.createCV("x", x),
                            datalogger.createCV("y", y),
                            datalogger.createCV("xya", xya),
                            datalogger.createCV("xyhi", xyhi))
                    }
                }
            }

            // projection in YZ plane...
            rsq = ysq + zsq
            yzRsq = rsq - Inertia * (rsq - yzRsq)
            yzlo = Math.min(yzlo, yzRsq)
            if (yzRsq > yzhi) {
                yzhi = yzRsq
                yza = Math.atan2(z, y) // angle (anticlockwise from Y axis)
                if ((stamp - yzLast) > MinPeakSeparation) { // need to clock new peak
                    yzPeaks.push(j)
                    yzLast = stamp
                    if (logging) {
                        datalogger.log(
                            datalogger.createCV("index", j),
                            datalogger.createCV("y", y),
                            datalogger.createCV("z", z),
                            datalogger.createCV("yza", yza),
                            datalogger.createCV("yzhi", yzhi))
                    }
                }
            }

            // projection in ZX plane...
            rsq = zsq + xsq
            zxRsq = rsq - Inertia * (rsq - zxRsq)
            zxlo = Math.min(zxlo, zxRsq)
            if (zxRsq > zxhi) {
                zxhi = zxRsq
                zxa = Math.atan2(x, z)  // angle (anticlockwise from Z axis)
                if ((stamp - zxLast) > MinPeakSeparation) { // need to clock new peak
                    zxPeaks.push(j)
                    zxLast = stamp
                    if (logging) {
                        datalogger.log(
                            datalogger.createCV("index", j),
                            datalogger.createCV("z", z),
                            datalogger.createCV("x", x),
                            datalogger.createCV("zxa", zxa),
                            datalogger.createCV("zxhi", zxhi))
                    }
                }
            }
        }
        // get actual radii
        xylo = Math.sqrt(xylo)
        yzlo = Math.sqrt(yzlo)
        zxlo = Math.sqrt(zxlo)
        xyhi = Math.sqrt(xyhi)
        yzhi = Math.sqrt(yzhi)
        zxhi = Math.sqrt(zxhi)

        // check average field-strength
        strength = Math.sqrt(strength / nSamples)
        if (strength < MarginalField) {
            return -2  // "FIELD STRENGTH TOO WEAK"
        }

        // compute eccentricities (the ratio of longest to shortest radii) from their squares
        // (defending against divide-by-zero errors if Spin-circle had projected exactly edge-on)
        let xye = xyhi / (xylo + 1)
        let yze = yzhi / (yzlo + 1)
        let zxe = zxhi / (zxlo + 1)

        // Select the "roundest" axis-pair to use (the one with lowest eccentricity) and adopt 
        // its normalisation offsets, tilt angle and eccentricity. Call the new axes U & V.
        if (xye < yze) { // not YZ
            if (xye < zxe) { // not ZX either: definitely use XY
                plane = "XY"
                uDim = Dim.X
                vDim = Dim.Y
                uOff = xOff
                vOff = yOff
                theta = xya
                scale = xye
                // the other (more eccentric) Ellipses are more reliable for deriving period
                period = setPeriod(yzPeaks,zxPeaks)
            }
        } else { // not XY: roundest is either YZ or ZX
            if (yze < zxe) {
                plane = "YZ"
                uDim = Dim.Y
                vDim = Dim.Z
                uOff = yOff
                vOff = zOff
                theta = yza
                scale = yze
                period = setPeriod(zxPeaks, xyPeaks)
            } else {
                plane = "ZX"
                uDim = Dim.Z
                vDim = Dim.X
                uOff = zOff
                vOff = xOff
                theta = zxa
                scale = zxe
                period = setPeriod(xyPeaks, yzPeaks)
            }
        }
        /**/
        basic.clearScreen()
        basic.pause(100)
        basic.showString(plane)
        basic.pause(300)

        // check that period was correctly found...
        if (period < 0) return period 
        
    
        // We have successfully set up the projection parameters. Now we need to relate them to North.
        // Take the average of seven new readings to get a stable fix on the current heading
        let uRaw = 0
        let vRaw = 0
        if (testing) { //arbitrarily choose 10th reading for North
            uRaw = scanData[10][uDim]
            vRaw = scanData[10][vDim]
        } else {
            // get a new position as the sum of seven readings
            uRaw = input.magneticForce(uDim)
            vRaw = input.magneticForce(vDim)
            for (let i = 0; i < 6; i++) {
                basic.pause(5)
                uRaw += input.magneticForce(uDim)
                vRaw += input.magneticForce(vDim)
            }
        }

        // record its projection angle on the Spin-Circle as the (global) fixed bias to North
        toNorth = project(uRaw, vRaw)


        /* no longer relevant (I think) but testing will prove... (delete later)
        // For a clockwise scan, the maths requires the U-dim to lead the V-dim by 90 degrees
        // From the point of view of a buggy spinning clockwise from ~NW, the North vector appears 
        // to rotate anticlockwise, passing the +V axis first, and then the -U axis.
        // Check the timings of their first limits and, if necessary, swap the major/minor dimensions:
         if (axes[uDim].time0 < axes[vDim].time0) {
            let temp = uDim
            uDim = vDim
            vDim = temp
        }

        // Also check the polarities of these first limits in case the microbit
        // is mounted backwards: we expect the first uVal<0 and the first vVal>0
        uFlip = -(axes[uDim].limit0 / Math.abs(axes[uDim].limit0)) // = -1 if uVal>0
        vFlip = axes[vDim].limit0 / Math.abs(axes[vDim].limit0)    // = -1 if vVal<0
        */
        if (logging) {
            datalogger.setColumnTitles("uDim", "vDim", "uOff", "vOff",
              "theta", "scale", "period", "toNorth", "strength")

            datalogger.log(
                datalogger.createCV("uDim", uDim),
                datalogger.createCV("vDim", vDim),
                datalogger.createCV("uOff", uOff),
                datalogger.createCV("vOff", vOff),
                datalogger.createCV("theta",theta),
                datalogger.createCV("scale",scale),
                datalogger.createCV("period", period),
                datalogger.createCV("toNorth", toNorth),
                datalogger.createCV("strength",strength))
        }


        // return average RPM of original scan    
        return 60000 / period
    }
  



    /**
     * Read the magnetometer and return the current heading of the buggy in degrees
     */
    //% block="heading" 
    //% inlineInputMode=inline 
    //% weight=70
    export function degrees(): number {
        // read the magnetometer (seven times) and return the current heading in degrees from North
        let uRaw = 0
        let vRaw = 0
        if (testing) {
            let index = 10 + (test * turn45) // (we chose tenth sample for North)
            uRaw = scanData[index][uDim]
            vRaw = scanData[index][vDim]
            test ++
            if (test > 7) test = 0 // roll round points of the compass
        } else {
            // get the new position as the sum of seven readings
            uRaw = input.magneticForce(uDim)
            vRaw = input.magneticForce(vDim)
            for (let i = 0; i < 6; i++) {
                basic.pause(5)
                uRaw += input.magneticForce(uDim)
                vRaw += input.magneticForce(vDim)
            }
        }

        // project reading from ellipse to circle,
        let onCircle = project(uRaw, vRaw)
        // relate it to North and convert to degrees
        let angle = 57.29578 * (onCircle - toNorth)
        // angle currently runs anticlockwise from U-axis: subtract it from 90 degrees 
        // to reflect through the diagonal, so making it run clockwise from the V-axis
        angle = 90 - angle
        // roll any negative values into the positive range [0...359]
        angle = (angle + 360) % 360

        return angle
    }


    //% block="set test mode: $turnOn"
    //% inlineInputMode=inline 
    //% weight=50
    export function setTestMode(turnOn: boolean) {
        testing = turnOn
    }

    //% block="set logging mode: $loggingOn"
    //% inlineInputMode=inline 
    //% weight=40
    export function setLogMode(loggingOn: boolean) {
        logging = loggingOn
    }

    export function isLogging(): boolean {
        return logging
    }



// UTILITY FUNCTIONS

    // use some sample data, while debugging...
    function simulateScan() {
        scanTimes = [51901, 51929, 51957, 51985, 52013, 52041, 52069, 52097, 52125, 52155, 52182, 52209, 52237, 52265, 52293, 52321, 52349, 52377, 52405, 52433, 52461, 52489, 52517, 52545, 52573, 52601, 52629, 52657, 52685, 52716, 52741, 52769, 52797, 52825, 52853, 52881, 52909, 52937, 52965, 52993, 53021, 53049, 53077, 53105, 53134, 53161, 53192, 53217, 53245, 53273, 53301, 53329, 53357, 53385, 53413, 53441, 53469, 53497, 53525, 53553, 53581, 53609, 53640, 53669, 53697, 53725, 53754, 53781, 53809, 53837, 53865, 53893, 53921, 53950, 53978, 54009, 54037, 54065, 54093, 54121, 54149, 54177, 54205, 54233, 54262, 54289, 54318, 54346, 54377, 54405, 54433, 54461, 54489, 54517, 54545, 54573, 54601, 54630, 54658, 54686, 54716, 54741, 54769, 54797, 54825, 54853, 54882, 54913, 54941, 54969, 54997, 55025, 55054, 55081, 55110, 55138, 55166, 55196, 55221, 55249, 55277, 55306, 55334, 55362, 55390, 55421, 55449, 55477, 55505, 55533, 55562, 55590, 55621, 55649, 55677, 55705, 55734, 55762, 55790, 55821, 55849, 55877, 55906, 55934, 55962, 55993, 56021, 56049, 56077, 56106, 56137, 56165, 56193, 56221, 56250, 56278, 56309, 56337, 56366, 56394, 56425, 56453, 56482, 56510, 56541, 56569, 56598, 56626, 56657, 56686, 56714, 56742, 56773, 56802, 56830, 56861, 56890, 56921, 56949, 56977, 57005, 57038, 57065, 57094, 57122, 57150, 57178, 57206, 57234, 57262, 57290, 57318, 57346, 57374, 57402, 57430, 57458, 57486, 57514, 57542, 57570, 57598, 57626, 57654, 57682, 57710, 57744, 57770, 57798, 57826, 57854, 57882, 57910, 57938, 57966, 57994, 58022, 58050, 58078, 58106, 58134, 58162, 58190, 58218, 58246, 58274, 58308, 58333, 58362, 58390, 58418, 58446, 58474, 58502, 58530, 58558, 58586, 58614, 58642, 58670, 58698, 58726, 58754, 58782, 58810, 58845, 58874, 58902, 58930, 58958, 58986, 59014, 59042, 59070, 59098, 59126, 59154, 59182, 59210, 59238, 59266, 59294, 59329, 59358, 59386, 59414, 59442, 59470, 59498, 59526, 59554, 59582, 59610, 59638, 59666, 59694, 59722, 59757, 59786, 59814, 59842, 59870, 59898, 59926, 59954, 59982, 60010, 60038, 60066, 60094, 60122, 60150, 60185, 60214, 60249, 60278, 60306, 60334, 60362, 60390, 60418, 60446, 60474, 60502, 60530, 60558, 60593, 60621, 60650, 60678, 60706, 60734, 60762, 60790, 60818, 60853, 60881, 60910, 60938, 60966, 60994, 61022, 61050, 61078, 61113, 61142, 61170, 61198, 61226, 61254, 61282, 61310, 61345, 61374, 61402, 61430, 61458, 61486, 61514, 61549, 61578, 61606, 61634, 61662, 61690, 61725, 61754, 61782, 61810, 61838, 61873, 61902, 61930, 61958, 61986, 62021, 62050, 62078, 62106, 62134, 62169, 62198, 62226, 62254, 62289, 62318, 62346, 62375, 62409, 62438, 62466, 62502, 62530, 62558, 62594, 62622, 62650, 62685, 62713, 62741, 62769, 62797, 62825, 62853, 62881, 62909, 62938, 62965, 62993, 63031, 63058, 63086, 63114, 63142, 63171, 63198, 63226, 63254, 63282, 63310, 63338, 63366, 63394, 63422, 63450, 63478, 63506, 63534, 63562, 63590, 63618, 63646, 63674, 63712, 63738, 63766, 63794, 63822, 63851, 63878, 63906, 63934, 63962, 63990, 64018, 64046, 64074, 64103, 64130, 64158, 64186, 64214, 64242, 64270, 64309, 64338, 64366, 64395, 64422, 64450, 64478, 64506, 64534, 64562, 64590, 64618, 64647, 64674, 64703, 64731, 64759, 64787, 64815, 64853, 64882, 64910, 64938, 64966, 64994, 65022, 65050, 65078, 65107, 65134, 65162, 65191, 65218, 65247, 65275, 65303, 65331, 65369, 65398, 65426, 65454, 65482, 65510, 65538, 65566, 65594, 65623, 65650, 65679, 65707, 65735, 65763, 65802, 65830, 65858, 65898, 65925, 65953, 65994, 66022, 66050, 66078, 66106, 66134, 66162, 66190, 66219, 66246, 66274, 66302, 66330, 66358, 66386, 66414, 66442, 66471, 66498, 66526, 66555, 66583, 66611, 66639, 66666, 66695, 66722, 66751]
        let xData = [-165.9, -163.5, -161.4, -157.5, -153.6, -150.45, -146.4, -143.4, -140.25, -136.35, -134.1, -130.95, -126.75, -124.35, -121.95, -118.8, -116.4, -113.55, -111.45, -110.55, -108.6, -105.9, -104.4, -102.9, -101.7, -100.8, -99.15, -98.1, -97.65, -97.5, -96.6, -94.95, -94.2, -94.35, -93.9, -93.3, -92.85, -92.85, -94.35, -94.95, -95.1, -95.7, -96.3, -97.2, -98.25, -98.4, -99, -100.2, -102.45, -105.75, -108.9, -111.15, -115.2, -119.25, -123.3, -126.9, -130.05, -133.95, -138.45, -142.8, -147, -151.35, -152.25, -153.45, -153.15, -154.35, -154.2, -153.3, -153.6, -156.3, -158.1, -162.75, -164.85, -168, -172.5, -175.35, -179.1, -183.9, -186.45, -190.35, -193.5, -196.2, -197.85, -199.65, -201, -202.95, -204.6, -205.05, -205.95, -209.1, -211.05, -212.25, -214.35, -215.55, -219.3, -222.6, -225.75, -229.5, -234.6, -238.5, -243.45, -247.5, -251.4, -254.7, -257.25, -260.4, -262.05, -265.05, -266.4, -267.45, -270.15, -272.55, -274.2, -276, -277.35, -280.05, -283.95, -285.15, -286.5, -287.55, -288.9, -289.8, -290.1, -289.65, -289.95, -290.7, -291.3, -292.95, -294.45, -295.65, -297, -297.9, -298.95, -300, -300.6, -301.35, -301.65, -302.4, -303.6, -305.4, -305.55, -306.9, -307.65, -308.85, -308.1, -307.95, -307.35, -308.7, -307.65, -306.45, -306.45, -307.35, -307.2, -306.9, -305.4, -305.85, -306.75, -306.9, -306.15, -305.85, -306.15, -306.45, -306, -304.95, -303.75, -303.9, -303.9, -302.85, -301.95, -300.75, -300.15, -300, -298.5, -297.9, -296.85, -296.1, -295.5, -294.75, -292.65, -291.75, -290.25, -287.85, -285.45, -282.9, -280.2, -277.95, -274.8, -270, -267.75, -265.2, -262.2, -259.65, -257.1, -254.25, -251.25, -247.95, -244.8, -242.85, -238.8, -235.5, -233.1, -231.9, -229.8, -227.55, -224.7, -223.5, -221.1, -218.55, -215.7, -213.6, -211.8, -208.95, -205.65, -202.5, -199.5, -197.4, -192.9, -189.45, -186.6, -183.45, -180, -176.55, -173.4, -170.7, -167.25, -164.25, -162, -160.2, -157.35, -153.45, -150.6, -148.2, -145.2, -141.45, -139.2, -136.2, -134.4, -132.45, -129.9, -129.3, -128.7, -126.15, -126.3, -124.65, -123.9, -123.45, -121.8, -120.75, -121.65, -121.05, -121.5, -121.8, -120.9, -119.7, -119.85, -117.15, -115.8, -114.9, -112.95, -111.45, -110.25, -107.25, -106.05, -104.25, -102, -100.5, -99.3, -98.25, -96.9, -95.7, -94.8, -93.9, -93.15, -92.7, -90.9, -90.45, -90.9, -90, -89.7, -89.1, -88.95, -88.65, -88.2, -86.4, -85.95, -85.35, -84.6, -83.85, -84.45, -84.3, -83.7, -84, -83.7, -82.8, -81.9, -80.85, -80.4, -81.45, -80.85, -79.8, -81.9, -83.55, -83.85, -85.35, -85.35, -85.95, -87, -86.25, -85.65, -85.95, -85.8, -86.85, -87.6, -89.25, -90.15, -90.6, -92.1, -93.45, -95.25, -96, -97.95, -99.6, -101.25, -104.1, -106.65, -108.9, -112.05, -113.85, -117.45, -120.9, -123.6, -126.15, -127.95, -131.7, -135.45, -137.55, -141.6, -145.35, -149.1, -153.9, -157.05, -160.2, -164.25, -167.55, -169.65, -171.3, -173.25, -175.2, -177.6, -180.3, -182.85, -185.7, -189.6, -193.65, -198.6, -204.3, -209.7, -214.8, -221.7, -229.2, -236.4, -243.45, -249, -255.3, -261.75, -266.1, -269.85, -272.7, -275.25, -277.05, -278.55, -280.05, -282, -282.9, -284.85, -285.6, -286.5, -287.1, -288.15, -290.1, -291.9, -292.2, -293.1, -295.8, -297.75, -299.25, -301.5, -301.5, -303.3, -304.8, -305.1, -305.4, -305.7, -304.05, -306.75, -307.2, -307.95, -307.95, -308.55, -308.1, -308.7, -306.45, -306.3, -306, -306.3, -307.2, -307.65, -307.8, -308.55, -308.4, -309.3, -309, -308.4, -307.95, -307.5, -306.45, -306.45, -304.8, -304.35, -303.15, -302.1, -300.9, -300.75, -298.95, -297, -294.15, -292.35, -290.7, -289.05, -287.7, -285.45, -283.5, -283.2, -280.8, -279, -277.5, -274.8, -273.45, -272.7, -271.35, -270.9, -269.1, -266.7, -264.45, -262.05, -257.55, -253.05, -249.3, -246.45, -244.05, -241.2, -238.35, -237.15, -235.8, -233.1, -229.5, -226.8, -224.85, -221.85, -218.4, -214.8, -211.2, -209.25, -204.6, -199.05, -195.15, -189.9, -184.2, -179.25, -173.1, -168.6, -164.25, -159.15, -155.7, -152.1, -148.35, -145.5, -142.05, -139.05, -136.5, -133.5, -132.3, -131.1, -128.55, -126.45, -124.95, -124.35, -124.05, -122.4, -120.75, -120.45, -120.3, -120.45, -119.55, -118.95, -117.6, -116.7, -115.5, -114, -112.05, -110.1, -108.45, -106.95, -105.45, -103.2, -102.3, -99.6, -97.8, -96, -94.5, -92.55, -90.45, -88.5, -88.95, -87.45, -86.55, -86.4, -86.55, -87.45, -87.15]
        let yData = [2.55, 2.55, 3.15, 4.5, 4.05, 4.35, 5.25, 4.8, 3.75, 4.2, 3.9, 4.2, 3.9, 3.6, 3.3, 3.9, 2.4, 2.25, 1.5, 0.6, -1.05, -2.85, -4.35, -5.4, -8.1, -8.7, -10.05, -10.65, -11.4, -13.05, -14.55, -15, -18.3, -20.7, -22.95, -25.2, -27.15, -30.45, -34.05, -36.6, -38.25, -40.8, -43.2, -47.1, -48.6, -51, -53.55, -57.15, -61.05, -64.8, -67.35, -71.85, -76.65, -81.45, -85.95, -90.45, -94.35, -98.25, -102.3, -105.9, -109.35, -112.05, -113.55, -115.2, -115.8, -115.2, -113.25, -112.2, -111.45, -111.75, -113.1, -114.45, -116.4, -120.3, -123.45, -126.75, -129.6, -131.55, -133.95, -135.6, -136.5, -136.95, -138.45, -139.5, -140.25, -141.9, -143.25, -144.15, -144.6, -144.6, -144.75, -145.65, -145.8, -146.1, -146.4, -147.9, -148.35, -148.8, -148.95, -148.8, -149.85, -150.75, -151.5, -152.4, -152.55, -152.4, -153.3, -153, -153.15, -151.65, -150.9, -152.25, -152.4, -152.25, -152.1, -151.35, -151.35, -151.2, -149.1, -148.2, -146.55, -145.8, -145.35, -146.1, -144.9, -144.9, -144.75, -144.45, -144, -142.95, -140.85, -141, -139.65, -139.35, -139.05, -137.7, -136.65, -136.35, -134.25, -133.95, -131.55, -130.05, -129.9, -129.3, -127.5, -127.05, -125.55, -125.25, -124.5, -123, -121.95, -120.75, -119.4, -119.25, -118.35, -117.75, -115.65, -113.85, -111.6, -110.7, -107.85, -106.35, -103.2, -101.55, -99.3, -98.25, -95.7, -94.35, -91.65, -91.05, -89.4, -89.1, -87.9, -87.45, -85.65, -85.2, -83.55, -82.2, -80.25, -78.15, -76.2, -74.4, -70.65, -68.1, -66, -61.95, -59.4, -54.9, -51.75, -50.1, -47.1, -44.1, -42.6, -41.1, -39.45, -37.2, -34.8, -32.7, -30.9, -28.65, -25.2, -23.4, -22.2, -21.45, -20.1, -18, -17.1, -17.4, -16.5, -15.3, -13.8, -13.05, -12.6, -10.8, -8.85, -7.65, -5.1, -2.85, -0.6, 0.75, 2.25, 3.6, 4.95, 5.7, 6.15, 7.05, 7.65, 8.25, 9, 9.6, 10.05, 11.4, 11.7, 13.05, 13.35, 15, 15.6, 15.75, 14.7, 14.55, 14.85, 14.7, 13.65, 13.65, 13.8, 14.4, 14.55, 14.1, 14.7, 14.55, 14.55, 15, 15, 14.85, 14.85, 14.1, 13.8, 13.5, 12.3, 11.55, 12.3, 12, 12.45, 12, 11.1, 10.8, 11.1, 10.65, 10.35, 10.65, 10.65, 10.95, 9.9, 9.3, 7.65, 6.6, 5.25, 4.2, 2.4, 1.95, 1.35, 0.9, 0.45, -1.5, -2.85, -3.45, -3.15, -4.95, -6.75, -8.7, -10.2, -11.4, -12.75, -15.6, -16.05, -16.65, -16.8, -17.55, -18.45, -19.65, -19.8, -21.45, -22.65, -25.05, -26.85, -28.35, -30.45, -32.4, -33.15, -35.55, -37.05, -38.25, -39.6, -39.6, -41.7, -45.3, -46.8, -48.9, -51.15, -54.15, -58.35, -61.05, -63.3, -66.3, -69, -71.7, -75, -77.85, -81, -84.6, -87.6, -91.5, -95.55, -97.8, -100.95, -103.8, -106.5, -108.75, -110.85, -112.65, -115.35, -117, -119.4, -121.65, -124.2, -126, -127.5, -129.45, -131.55, -132.9, -133.8, -135, -136.95, -139.2, -141.3, -142.8, -144.45, -147.15, -149.85, -151.05, -152.25, -152.85, -153.3, -154.8, -154.05, -154.05, -154.8, -154.65, -153.75, -153.15, -152.4, -152.55, -151.95, -150.6, -150.3, -149.55, -149.7, -149.4, -148.2, -147.3, -147.3, -146.55, -146.25, -144.75, -143.1, -142.05, -140.55, -138, -135.9, -135, -133.95, -133.05, -131.1, -130.8, -130.35, -129.9, -129.3, -128.4, -127.35, -127.35, -125.7, -124.65, -123.75, -122.55, -122.7, -121.65, -120.75, -119.7, -119.1, -118.35, -116.85, -114.15, -112.95, -111.3, -110.4, -108.75, -106.8, -105, -103.65, -101.4, -98.4, -94.8, -91.05, -87.9, -84.6, -81.75, -78.3, -76.35, -74.7, -73.8, -71.7, -69.75, -68.55, -67.35, -65.85, -63.6, -62.7, -61.95, -61.2, -59.1, -57.9, -56.4, -54.45, -51.3, -48.6, -45.3, -42.6, -39.3, -35.85, -34.35, -32.55, -31.65, -30.75, -28.65, -28.2, -27.3, -25.65, -23.55, -20.85, -19.8, -18.45, -15.75, -14.4, -12, -9, -6.15, -2.4, -0.15, 1.65, 3.75, 5.55, 5.55, 6.45, 7.05, 8.4, 9.75, 10.95, 11.7, 13.8, 14.25, 14.55, 13.95, 14.25, 14.1, 13.95, 13.5, 14.4, 14.1, 14.4, 13.35, 13.35, 13.05, 12.15, 10.95, 11.25, 11.7, 11.7, 12.15, 11.85, 12.75, 12.9, 13.05, 12.9, 12.45, 11.4, 10.8, 10.5, 8.85, 7.5, 6.15, 4.95, 4.35, 3.15, 1.05, 0.6, -1.05, -2.7, -3, -5.25, -6.3]
        let zData = [379.2, 377.7, 374.7, 369.15, 366.15, 362.85, 359.25, 356.85, 353.4, 349.8, 346.8, 343.65, 338.85, 334.8, 330.6, 325.5, 321, 316.8, 312.9, 310.2, 307.95, 303, 300.9, 297.45, 294.9, 291.6, 289.65, 286.65, 285.15, 282.3, 280.8, 277.95, 274.8, 269.85, 265.5, 262.5, 258.6, 254.4, 252.45, 249.6, 247.95, 247.05, 245.85, 244.65, 243.9, 242.55, 242.1, 241.35, 240.45, 239.55, 237.45, 236.55, 234.6, 234, 233.7, 233.85, 233.7, 235.65, 236.1, 238.35, 238.35, 240.15, 240.75, 240.6, 240.75, 241.65, 240.9, 241.2, 241.2, 241.8, 244.2, 245.55, 247.65, 249.9, 253.65, 255.15, 258.45, 260.25, 262.65, 263.7, 265.8, 266.25, 267.9, 267.9, 269.1, 270.45, 272.55, 273, 274.95, 275.7, 277.8, 279.3, 281.4, 283.65, 287.4, 290.1, 294.9, 298.2, 302.25, 305.85, 309, 312.6, 316.2, 319.35, 322.65, 326.4, 328.8, 332.25, 333.9, 337.05, 338.25, 339.75, 340.35, 342, 343.5, 345.75, 346.95, 350.1, 352.8, 356.1, 358.05, 359.55, 362.25, 364.65, 366.3, 368.7, 369.45, 371.85, 373.8, 374.85, 375.75, 377.85, 378.45, 380.25, 382.8, 385.35, 387.45, 388.65, 390.3, 391.2, 392.4, 392.55, 393.75, 394.65, 397.5, 397.95, 400.95, 402.15, 403.8, 403.5, 405.3, 406.35, 408.3, 409.2, 410.85, 412.05, 415.35, 416.85, 418.05, 419.1, 421.5, 423.75, 426.3, 427.05, 428.1, 429.75, 431.4, 431.1, 432.3, 431.4, 433.05, 433.8, 435.15, 435, 435.3, 435.45, 437.25, 437.1, 437.1, 436.35, 436.8, 437.7, 437.7, 437.7, 437.85, 438.3, 438, 438.15, 437.4, 436.65, 436.05, 435.3, 433.65, 433.65, 432.6, 432.15, 431.25, 430.35, 429.6, 429.6, 427.05, 425.85, 424.35, 422.7, 421.5, 420.6, 418.65, 418.8, 417.45, 416.55, 416.55, 414.45, 412.8, 411.3, 408.3, 406.35, 403.8, 400.95, 398.55, 396, 393.6, 390.6, 388.8, 386.7, 383.7, 381.6, 378.3, 376.2, 374.85, 371.4, 368.55, 365.7, 363, 360.6, 357.75, 355.5, 353.25, 351.75, 350.55, 349.05, 348.75, 347.25, 345.9, 345.45, 345, 345.15, 344.1, 342.15, 342.45, 342.15, 342, 340.35, 338.55, 337.8, 337.2, 335.25, 333.75, 331.65, 329.25, 327, 324.75, 321.75, 319.35, 317.55, 314.85, 312.6, 310.5, 308.4, 306.6, 305.25, 302.1, 301.5, 300.9, 299.55, 298.35, 297.3, 296.1, 294.9, 291.45, 288.9, 286.5, 285.15, 283.8, 280.65, 278.7, 277.5, 276.6, 274.35, 272.7, 269.85, 268.8, 267.45, 265.65, 263.55, 262.8, 260.7, 260.85, 259.65, 258.15, 257.4, 256.35, 254.25, 252.45, 249.75, 249.3, 248.4, 247.5, 246.45, 245.85, 244.8, 244.2, 241.95, 240.6, 239.7, 238.2, 237.9, 237, 235.5, 234.45, 233.25, 232.2, 232.05, 231.45, 231.6, 232.2, 231.45, 231.75, 231.45, 229.95, 230.4, 228.9, 228.9, 231, 232.35, 232.2, 234.75, 235.65, 238.5, 238.95, 240.3, 241.35, 244.05, 246.45, 248.4, 249.6, 252.3, 252.6, 253.65, 254.85, 254.7, 256.35, 259.05, 261.6, 266.1, 269.85, 274.95, 279.9, 284.7, 289.5, 295.35, 299.7, 305.4, 310.2, 316.05, 322.5, 327.75, 333, 337.65, 340.8, 344.7, 347.4, 349.05, 351.45, 352.95, 354.3, 356.7, 358.2, 360.75, 361.95, 363.9, 365.55, 368.25, 370.05, 372.75, 375.45, 379.05, 382.5, 385.05, 387.3, 390.45, 393.15, 394.35, 396.45, 396.9, 398.4, 400.35, 401.85, 401.85, 403.65, 403.65, 404.1, 404.7, 404.55, 405.45, 406.5, 406.95, 407.55, 409.35, 410.55, 412.8, 413.25, 414.15, 416.1, 416.7, 417.45, 418.95, 420, 421.2, 423.3, 424.35, 427.8, 430.35, 431.7, 434.25, 436.8, 438.45, 439.2, 440.1, 440.4, 441.45, 440.25, 439.95, 439.5, 440.1, 440.25, 440.7, 440.85, 442.35, 442.35, 442.35, 442.5, 441.6, 441.3, 439.8, 438.6, 437.4, 436.2, 434.1, 433.05, 430.5, 429.9, 428.25, 427.05, 426.9, 425.7, 423.9, 423.45, 421.05, 418.5, 416.1, 413.7, 410.85, 408.45, 405.6, 403.05, 400.65, 398.25, 394.65, 391.8, 388.5, 384.6, 381.6, 378.6, 375.45, 371.1, 368.4, 366, 363.6, 361.35, 359.25, 356.7, 354.6, 353.1, 350.85, 350.1, 348.75, 347.4, 346.5, 347.1, 345.45, 344.55, 343.05, 340.65, 339, 337.2, 335.7, 334.05, 332.25, 328.2, 325.95, 323.55, 321.6, 318, 314.55, 311.85, 310.2, 307.95, 305.55, 302.25, 300.3, 298.5, 295.8, 294.3, 291.9, 288.6, 286.05, 283.2, 281.1, 277.35, 274.35]
        
        // transpose three arrays into array of triples
        for (let n = 0; n < scanTimes.length; n++){
            let xyz = []
            xyz.push(xData[n])
            xyz.push(yData[n])
            xyz.push(zData[n])
            scanData.push(xyz)
        }
    }

    // Perform periodicity analysis from two sets of peaks....
    function setPeriod(peaks1: number[], peaks2: number[]): number {
        // did we spin enough to give at least three peak values in both sets?
        if ((peaks1.length < 3) || (peaks2.length < 3))  {
            return -3 // NOT ENOUGH SCAN ROTATION
        }
        // An Ellipse has two peaks, so period is twice the average time between maximae
        let first = peaks1[0]
        let last = peaks1.pop()
        let gaps = peaks1.length //...having popped the last one
        // split time between gaps
        let period1 = (scanTimes[last] - scanTimes[first]) / gaps

        first = peaks2[0]
        last = peaks2.pop()
        gaps = peaks2.length
        let period2 = (scanTimes[last] - scanTimes[first]) / gaps

        period = period1 + period2 // effectively averaging the two results

        /************** global just for testing purposes *************/
        turn45 = Math.floor((last - first) / (4 * gaps)) // ~ #samples covering an octant
        /*************************************************************/

        return period
    }


    // Transform a point on the off-centre projected Ellipse back onto the centred Spin-Circle 
    // and return its angle (in radians) anticlockwise from the horizontal U-axis
    // Uses global projection parameters: {uDim, vDim, uOff, vOff, theta, scale}
    function project(uRaw: number, vRaw: number): number {
        // shift to start the vector at the origin
        let u = uRaw - uOff
        let v = vRaw - vOff
        // rotate by the major-axis angle theta (check direction!)
        let uNew = u * Math.cos(theta) - v * Math.sin(theta)
        let vNew = u * Math.sin(theta) + v * Math.cos(theta)
        // scale up V to match U
        let vScaled = vNew * scale
        // return projected angle (undoing the rotation we just applied)
        let angle = Math.atan2(vScaled, uNew) - theta
        if (logging) {

            datalogger.setColumnTitles("uRaw", "vRaw", "u", "v", "uNew", "vNew", "vScaled", "angle")
            datalogger.log(datalogger.createCV("uRaw", uRaw),
                datalogger.createCV("vRaw", vRaw),
                datalogger.createCV("u", u),
                datalogger.createCV("v", v),
                datalogger.createCV("uNew", uNew),
                datalogger.createCV("vNew", vNew),
                datalogger.createCV("vScaled", vScaled),
                datalogger.createCV("angle", angle) )
        }
        
        return angle
    }
    
}
