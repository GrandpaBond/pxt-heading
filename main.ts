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

    let scanTimes: number[] = [] // sequence of time-stamps for scanned readings 
    let scanData: number[][]= [] // scanned sequence of [X,Y,Z] magnetometer readings  
    let uDim = -1 // the "horizontal" axis (pointing East) for transformed readings (called U)
    let vDim = -1 // the "vertical" axis (pointing North) for transformed readings (called V)
    let uOff = 0 // the offset needed to re-centre readings along the U-axis
    let vOff = 0 // the offset needed to re-centre readings along the V-axis
    let theta = 0 // the angle to rotate readings so the projected ellipse aligns with the U & V axes
    let scale = 0 // the scaling to stretch V-axis readings from the ellipse onto a circle
    let toNorth = 0 // the angular bias to be added (so that North = 0)
    let strength = 0 // the average magnetic field-strength observed by the magnetometer
    let uFlip = 1 // set to -1 if uDim polarity is inverted
    let vFlip = 1 // set to -1 if vDim polarity is inverted
    let plane = "??" // the projection plane we are using: "XY","YZ" or "ZX"
    let testing = false  // test mode flag
    let test = 0         // selector for test sample

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
            return
        }

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
    /* 
    As the buggy spins, the magnetic field-vector sweeps out a Spin-Circle on the surface of a sphere.
    In the fully general case, this projects onto the plane of each pair of orthogonal axes (XY,YZ,ZX) 
    as an ellipse with a certain eccentricity. We will get the best heading discrimination from the 
    plane with the least eccentric ellipse, and having selected those two axes, we'll need to 
    transform readings around the ellipse so that they lie on a circle, giving a relative 
    angle that can (eventually) be offset by a fixed bias to return the true heading.
    */

        // we need at least 2 second's worth of scanned readings...
        if (scanTimes.length < 80) {
            return -1 // "NOT ENOUGH SCAN DATA"
        }
        
        // Find the raw extrema, and note when they happened by detecting change of sign of slope.
        // Data are already rolling sums of 7 values, so subtract 4 behind from 4 ahead to get slope 
        // Record in peaks[][] the index of each inflection point, to work out rotation rate later
        // (but ignore repeated inflections if they are too close in time)    
        let xlo = 999
        let ylo = 999
        let zlo = 999
        let xhi = -999
        let yhi = -999
        let zhi = -999
        // peaks is an array [X,Y,Z] of three sub-arrays of axis maximae, added as we discover them
        //(NOTE: the sub-array lengths will often differ between dimensions)
        let peaks: number[][] = [[],[],[]]
        // get the slopes for the fourth sample
        let dx = scanData[7][Dim.X] - scanData[0][Dim.X]
        let dy = scanData[7][Dim.Y] - scanData[0][Dim.Y]
        let dz = scanData[7][Dim.Z] - scanData[0][Dim.Z]
        let xLast = scanTimes[0]
        let yLast = scanTimes[0]
        let zLast = scanTimes[0]
        let dxWas = 0
        let dyWas = 0
        let dzWas = 0
        let v = 0
        // loop over central samples, (skipping 4-sample "wings" at each end)
        for (let i = 4; i < scanTimes.length - 5; i++) {
            dxWas = dx
            dx = scanData[i + 4][Dim.X] - scanData[i - 4][Dim.X]
            v = scanData[i][Dim.X]
            if ((dx * dxWas) < 0) // scanData[i] is a local max or min for X
            {
                if (dxWas > 0) {
                    xhi = v // reached local maximum
                    if (scanTimes[i] - xLast > MinPeakSeparation) {
                        peaks[Dim.X].push(i)
                        xLast = scanTimes[i]
                    }
                } else xlo = v // reached local minimum
            }

            dyWas = dy
            dy = scanData[i + 4][Dim.Y] - scanData[i - 4][Dim.Y]
            v = scanData[i][Dim.Y]
            if ((dy * dyWas) < 0) // scanData[i] is a local max or min for Y
            {
                if (dyWas > 0) {
                    yhi = v // reached local maximum
                    if (scanTimes[i] - yLast > MinPeakSeparation) {
                        peaks[Dim.Y].push(i)
                        yLast = scanTimes[i]
                    }
                } else ylo = v // reached local minimum
            }

            dzWas = dz
            dz = scanData[i + 4][Dim.Z] - scanData[i - 4][Dim.Z]
            v = scanData[i][Dim.Z]
            if ((dz * dzWas) < 0) // scanData[i] is a local max or min for Z
            {
                if (dzWas > 0) {
                    zhi = v // reached local maximum
                    if (scanTimes[i] - zLast > MinPeakSeparation) {
                        peaks[Dim.Z].push(i)
                        zLast = scanTimes[i]
                    }
                } else zlo = v // reached local minimum
            }
        }

     // just use the latest extremes to set the normalisation offsets
        let xOff = (xhi + xlo) / 2
        let yOff = (yhi + ylo) / 2
        let zOff = (zhi + zlo) / 2

        // now find the extreme radii-squared for each axis-pair 
        let xylo = 99999
        let yzlo = 99999
        let zxlo = 99999
        let xyhi = -99999
        let yzhi = -99999
        let zxhi = -99999
        let x = 0
        let y = 0
        let z = 0
        let xsq = 0
        let ysq = 0
        let zsq = 0
        let rsq = 0
        let xya = 0
        let yza = 0
        let zxa = 0
        let nSamples = scanTimes.length

        for (let j = 0; j < nSamples; j++) {
            // extract next readings and normalise them
            x = scanData[j][Dim.X] - xOff
            y = scanData[j][Dim.Y] - yOff
            z = scanData[j][Dim.Z] - zOff
            xsq = x * x
            ysq = y * y
            zsq = z * z
            strength += xsq + ysq + zsq // accumulate square of field strengths

            // projection in XY plane
            rsq = xsq + ysq // Pythagoras: radius-squared = sum-of-squares
            if (rsq < xylo) xylo = rsq // shortest so far...
            if (rsq > xyhi) {
                xyhi = rsq // longest so far...
                xya = Math.atan2(x, y) // angle (anticlockwise from X axis)
            }

            // projection in YZ plane
            rsq = ysq + zsq
            if (rsq < yzlo) yzlo = rsq
            if (rsq > yzhi) {
                yzhi = rsq
                yza = Math.atan2(y, z) // angle (anticlockwise from Y axis)
            }

            // projection in ZX plane
            rsq = zsq + xsq
            if (rsq < zxlo) zxlo = rsq
            if (rsq > zxhi) {
                zxhi = rsq
                zxa = Math.atan2(z, x)  // angle (anticlockwise from Z axis)
            }
        }

        // check average field-strength
        strength = Math.sqrt(strength / nSamples)
        if (strength < MarginalField) {
            return -2  // "FIELD STRENGTH TOO WEAK"
        }

        // compute eccentricities from latest max & min radii-squared
        // (defending against divide-by-zero errors if Spin-circle was exactly edge-on)
        let xye = Math.sqrt(xyhi / (xylo + 0.0001))
        let yze = Math.sqrt(yzhi / (yzlo + 0.0001))
        let zxe = Math.sqrt(zxhi / (zxlo + 0.0001))
    
        // select best axis-pair to use (the one with lowest eccentricity)
        if (xye < yze) { // not YZ
            if (xye < zxe) { // not ZX either, so use XY
                plane = "XY"
                uDim = Dim.X
                vDim = Dim.Y
                uOff = xOff
                vOff = yOff
                theta = xya
                scale = xye
            }
        } else { // not XY: smallest is either YZ or ZX
            if (yze < zxe) {
                plane = "YZ"
                uDim = Dim.Y
                vDim = Dim.Z
                uOff = yOff
                vOff = zOff
                theta = yza
                scale = yze
            } else {
                plane = "ZX"
                uDim = Dim.Z
                vDim = Dim.X
                uOff = zOff
                vOff = xOff
                theta = zxa
                scale = zxe
            }
        }

        // did we spin enough to give at least a couple of peak values?
        if (peaks[uDim].length < 2) {
            return -3 // NOT ENOUGH SCAN ROTATION
        }

        // period is average time between U-axis maximae
        let first = peaks[uDim][0]
        let last = peaks[uDim].pop()
        let gaps = peaks[uDim].length //... now that we've popped the last one
        // split time between gaps
        let period = (scanTimes[last] - scanTimes[first]) / gaps


// We have successfully set up the projection parameters. Now we need to relate them to North.
// Take the average of seven new readings to get a stable fix on the current heading
        let uRaw = 0
        let vRaw = 0
        if (testing) { // choose some arbitrary reading for North (using X for U; Z for V)
            uRaw = scanData[40][Dim.X]
            vRaw = scanData[40][Dim.Z]
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

        // record its projection angle as the (global) fixed bias to North
        toNorth = project(uRaw, vRaw)


 /* no longer relevant (I think) but check sometime --delete later
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
        // set up datalogger for subsequent calls to project() from heading.degrees()
        datalogger.deleteLog()
        datalogger.includeTimestamp(FlashLogTimeStampFormat.None)
        datalogger.setColumnTitles("uRaw", "vRaw", "u", "v", "uNew", "vNew", "vScaled", "angle")

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
        if (testing) { // NOTE: a-priori knowledge: U=X & V=Z for test data!
            uRaw = scanData[test][Dim.X]
            vRaw = scanData[test][Dim.Z]
            test += 4
            if (test > scanTimes.length - 5) test = 0 // roll round before running out
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

        // project reading from ellipse to circle, relating it to North and converting to degrees
        let angle = 57.29578 * (project(uRaw, vRaw) - toNorth)
        // angle currently runs anticlockwise from U-axis: subtract it from 90 degrees 
        // to reflect through the diagonal, so making it run clockwise from the V-axis
        angle = 90 - angle
        // roll any negative values into the positive range [0...359]
        angle = (angle + 360) % 360

        return angle
    }


    //% block="switch test mode" 
    //% inlineInputMode=inline 
    //% weight=50
    export function testMode(turnOn: boolean) {
        testing = turnOn
    } 

// UTILITY FUNCTIONS

    // use some sample data, while debugging...
    function simulateScan() {
        scanTimes = [23673, 23701, 23729, 23757, 23785, 23813, 23841, 23869, 23899, 23925, 23953, 23981, 24009, 24037, 24065, 24093, 24121, 24149, 24177, 24205, 24233, 24261, 24290, 24317, 24345, 24373, 24401, 24429, 24460, 24485, 24513, 24541, 24569, 24597, 24626, 24653, 24681, 24709, 24737, 24765, 24793, 24821, 24849, 24877, 24905, 24936, 24961, 24989, 25017, 25045, 25073, 25101, 25129, 25157, 25185, 25213, 25241, 25270, 25298, 25326, 25353, 25384, 25409, 25437, 25465, 25493, 25521, 25549, 25577, 25605, 25634, 25662, 25690, 25718, 25748, 25773, 25801, 25829, 25857, 25885, 25913, 25941, 25970, 25998, 26026, 26054, 26082, 26112, 26137, 26165, 26193, 26221, 26250, 26278, 26306, 26334, 26362, 26390, 26421, 26449, 26477, 26505, 26533, 26562, 26590, 26618, 26649, 26677, 26705, 26733, 26761, 26790, 26818, 26846, 26874, 26902, 26932, 26958, 26985, 27014, 27042, 27070, 27098, 27126, 27157, 27185, 27214, 27241, 27270, 27298, 27326, 27357, 27385, 27413, 27441, 27470]
        let xData = [-32.4, -23.55, -8.85, 4.5, 22.35, 38.4, 55.2, 71.4, 85.95, 99.15, 112.65, 126, 138, 148.2, 157.35, 164.1, 168, 170.1, 169.65, 165.45, 158.1, 150.75, 142.8, 134.1, 123.45, 110.4, 95.85, 81, 65.1, 47.85, 29.7, 10.95, -5.85, -21.3, -34.95, -48.75, -61.65, -71.85, -82.35, -93.3, -101.7, -109.05, -114.9, -117.45, -120.75, -118.2, -114.9, -109.5, -101.1, -92.4, -82.95, -70.35, -60.6, -46.2, -29.7, -13.2, 4.05, 23.85, 41.25, 61.05, 80.85, 97.5, 113.1, 129.45, 140.85, 151.05, 159.9, 164.7, 168, 168.45, 165.6, 164.25, 161.55, 156.15, 149.4, 141, 132.3, 121.2, 107.85, 91.2, 74.7, 56.55, 39.3, 20.4, 2.1, -14.7, -28.65, -43.05, -54.3, -66.6, -77.55, -86.55, -96.75, -105.6, -112.05, -118.8, -122.4, -123.6, -123.6, -120.15, -112.95, -105, -93, -81.6, -67.8, -52.95, -37.65, -21, -3.15, 15, 34.2, 52.8, 70.35, 87.6, 102.9, 117, 128.25, 139.05, 148.05, 156.3, 163.35, 167.85, 170.4, 171.6, 170.85, 166.65, 160.65, 153.6, 143.1, 131.25, 115.5, 100.5, 81.45, 63.45, 43.2, 24.6]
        let yData = [-853.2, -849.9, -850.65, -850.5, -850.65, -850.65, -849.6, -849.9, -853.35, -852, -851.55, -851.4, -851.55, -851.85, -851.1, -848.85, -849.45, -851.25, -851.55, -852.3, -851.85, -854.1, -856.5, -858.75, -859.05, -859.2, -857.55, -861.15, -860.1, -859.2, -857.25, -857.1, -857.85, -859.95, -858.3, -857.7, -857.7, -858.15, -858, -856.2, -856.65, -856.05, -856.8, -855.75, -855.6, -855.45, -855.9, -854.1, -853.05, -851.4, -849.75, -849.6, -848.25, -847.5, -847.5, -846.6, -846.75, -848.55, -849, -849.45, -848.25, -847.5, -847.05, -847.05, -844.95, -845.55, -844.95, -848.55, -850.65, -851.25, -851.25, -854.55, -852.9, -853.95, -854.7, -855, -856.65, -858.15, -857.7, -859.05, -859.95, -858.6, -859.65, -860.55, -862.95, -862.95, -862.35, -862.95, -863.55, -861.9, -860.85, -857.55, -857.1, -856.8, -856.2, -853.95, -851.85, -850.5, -850.95, -849.3, -849.45, -848.25, -849.45, -849.6, -850.05, -849.3, -848.4, -847.5, -847.5, -846.9, -846.75, -846.3, -847.35, -847.5, -846.45, -847.65, -848.1, -848.7, -847.95, -846.9, -847.8, -849.6, -848.85, -848.55, -850.05, -852.75, -853.2, -855, -857.4, -856.35, -856.35, -856.8, -857.7, -857.4, -856.95, -857.55]
        let zData = [-195.75, -193.05, -189.45, -188.1, -185.55, -183.9, -186.6, -191.7, -196.5, -204.45, -214.05, -223.8, -236.85, -251.1, -265.8, -283.8, -298.8, -314.55, -333.45, -349.35, -366, -382.65, -397.35, -412.8, -425.55, -438, -448.5, -454.5, -458.55, -459.9, -462.15, -461.7, -457.95, -452.4, -446.55, -438.75, -430.35, -419.4, -405.6, -393.3, -379.2, -364.05, -349.8, -335.55, -319.05, -303.15, -288.6, -274.05, -258.6, -244.65, -233.25, -221.85, -215.1, -204.6, -196.65, -192.15, -187.35, -184.05, -185.7, -186, -190.5, -197.55, -208.05, -219, -230.85, -244.35, -259.5, -276, -293.1, -309.45, -328.2, -344.7, -361.65, -375.6, -390.9, -406.2, -419.1, -429.75, -441, -447.3, -457.2, -462.6, -464.1, -465.9, -463.2, -460.05, -453.9, -446.85, -436.8, -427.2, -414.9, -401.55, -387.9, -376.05, -360, -347.55, -331.65, -316.05, -301.5, -286.65, -270.15, -256.2, -238.8, -226.95, -213.9, -204, -195.75, -190.65, -185.4, -182.55, -181.8, -183.15, -186.3, -190.8, -197.55, -207, -218.25, -232.05, -246.15, -263.1, -279.15, -295.2, -313.5, -332.4, -348.15, -366.45, -382.8, -398.7, -413.1, -426.15, -438.6, -447.15, -453.9, -459.9, -462.75, -462.9]
   // transpose three arrays into array of triples
        for (let n = 0; n < scanTimes.length; n++){
            let xyz = []
            xyz.push(xData[n])
            xyz.push(yData[n])
            xyz.push(zData[n])
            scanData.push(xyz)
        }
    }

    export function dumpData() {
        datalogger.deleteLog()
        datalogger.includeTimestamp(FlashLogTimeStampFormat.None)
        datalogger.setColumnTitles("t", "x", "y", "z")
        for (let i = 0; i < scanTimes.length; i++) {
            datalogger.log(datalogger.createCV("t", scanTimes[i]),
                datalogger.createCV("x", scanData[i][Dim.X]),
                datalogger.createCV("y", scanData[i][Dim.Y]),
                datalogger.createCV("z", scanData[i][Dim.Z]))
        }
    }


    // Transform a point on the off-centre projected ellipse back onto a centred circle of headings 
    // and return its angle (in radians) anticlockwise from the horizontal U-axis
    function project(uRaw: number, vRaw: number): number {
        // shift to start the vector at the origin
        let u = uRaw - uOff
        let v = vRaw - vOff
        // rotate by the major-axis angle theta (check direction!)
        let uNew = u * Math.cos(theta) - v * Math.sin(theta)
        let vNew = u * Math.sin(theta) + v * Math.cos(theta)
        // scale up V to match U
        let vScaled = vNew * scale
        // return projected angle
        let angle = Math.atan2(vScaled, uNew)

        datalogger.log(datalogger.createCV("uRaw", uRaw),
            datalogger.createCV("vRaw", vRaw),
            datalogger.createCV("u", u),
            datalogger.createCV("v", v),
            datalogger.createCV("uNew", uNew),
            datalogger.createCV("vNew", vNew),
            datalogger.createCV("vScaled", vScaled),
            datalogger.createCV("angle", angle) )
        
        return angle
    }
    
}
