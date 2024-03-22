/**
 * An extension providing a compass-heading for a buggy located anywhere on the globe 
 * (except at the magnetic poles!), with any mounting orientation for the microbit.
 */

//% color=#6080e0 weight=40 icon="\uf14e" block="Heading" 
namespace heading {

    // ENUMERATIONS

    enum Dim { // ...just for brevity!
        X = Dimension.X,
        Y = Dimension.Y,
        Z = Dimension.Z
    }

    enum View { // the axis-pairs spanning the three possible projection planes
        XY,
        YZ,
        ZX
    }

    // CONSTANTS

    const MarginalField = 30 // minimum acceptable field-strength for magnetometer readings
    const Inertia = 0.95 // ratio of old to new radius readings (for inertial smoothing)
    const Window = 100 // minimum spacing (in samples) of Ellipse major-axis candidates

    // CLASSES

    // Characteristics of the Ellipse formed when projecting the Spin-Circle onto a View plane
    class Ellipse {
        plane: string; // name (for debug)
        uDim: number; // horizontal axis of this View
        vDim: number; // vertical axis of this View
        uOff: number; // horizontal offset needed to re-centre the Ellipse along the U-axis
        vOff: number; // vertical offset needed to re-centre the Ellipse along the V-axis
        // properties derived from the Spin-Circle scan:
        rSq: number;  // smoothed radius-squared
        rSqWas: number; // previous smoothed radius-squared
        hiRsq: number; // max radius-squared 
        loRsq: number; // min radius-squared
        peaks: number[]; // array of indices of maxima
        // rotation params for aligning Ellipse major axis horizontally:
        theta: number; // angle (in radians) anticlockwise from U-axis to major axis 
        cosTheta: number;
        sinTheta: number;
        scale: number; // stretch in V needed to make Ellipse circular again 
        semiPeriod: number; // this View's assessment of half of the rotation time

        constructor(plane: string, dim0: number, dim1: number, off0: number, off1: number) {
            this.plane = plane // (as a DEBUG aid)
            this.uDim = dim0
            this.vDim = dim1
            this.uOff = off0
            this.vOff = off1
            this.peaks = []       
        }

        // "The square on the hypotenuse..."
        pythagoras(uRaw: number, vRaw: number): number {
            let u = uRaw - this.uOff
            let v = vRaw - this.vOff
            let r = (u * u) + (v * v)
            return r
        }


        // initialise the smoothing history
        firstSample(uRaw: number, vRaw:number) {
            this.rSq = this.pythagoras(uRaw,vRaw)
            this.rSqWas = this.rSq
            this.hiRsq = this.rSq
            this.loRsq = this.rSq
        }

        // process another scan sample
        nextSample(index: number, uRaw: number, vRaw: number) {
            let previous = -Window // always permit first peak
            // while tracking the radius, we use inertial smoothing to reduce multiple
            // spurious peak-detections due to minor fluctuations in readings
            let smooth = (this.rSq * Inertia) + (this.pythagoras(uRaw, vRaw) * (1 - Inertia))
            this.hiRsq = Math.max(this.hiRsq, smooth) // longest so far...
            this.loRsq = Math.min(this.loRsq, smooth) // shortest so far...
            // need to clock new peak if slope changes from rising to falling, but not too recently
            if (   (this.rSq >= this.rSqWas) 
                && (this.rSq >= smooth)
                && ((index - previous) >= Window) )  {
                this.peaks.push(index) // (technically, previous sample was the peak, but we're near enough)
                if (logging) {
                    datalogger.log(
                        datalogger.createCV("view", this.plane),
                        datalogger.createCV("index", index),
                        datalogger.createCV("loRsq", this.loRsq),
                        datalogger.createCV("hiRsq", this.hiRsq))

                }
                previous = index
            }
            this.rSqWas = this.rSq
            this.rSq = smooth
        }

        // analyse detected peaks[] to set transform rotation and scale; and measure period too
        analyse() {
            let hiR = Math.sqrt(this.hiRsq)
            let loR = Math.sqrt(this.loRsq)
            this.scale = hiR / (loR + 0.001) // == the eccentricity of this View's Ellipse
            // (defend against divide-by-zero errors, if Spin-circle had projected exactly edge-on)
            // Find angle of Ellipse's major-axis (anticlockwise from U-Axis) 
            let last = this.peaks.pop() // sample index of latest candidate
            let u = scanData[last][this.uDim] - this.uOff
            let v = scanData[last][this.vDim] - this.vOff
            this.theta = Math.atan2(u, v) //
            // we might as well remember these...
            this.cosTheta = u / hiR
            this.sinTheta = v / hiR
            // While we're at it, analyse the rotation period, in case anyone's interested
            // Ellipse's major-axis peaks fall 180 degrees apart...
            let first = this.peaks[0]
            let gaps = this.peaks.length //(having popped the [last] one)
            // ...so average time between peaks is for just half a rotation
            this.semiPeriod = (scanTimes[last] - scanTimes[first]) / gaps

            if (logging) {
                // prepare for analysis
                datalogger.setColumnTitles("view", "hiR", "loR", "first", "last", "semiPeriod")
                datalogger.log(
                    datalogger.createCV("view", this.plane),
                    datalogger.createCV("hiR", hiR),
                    datalogger.createCV("loR", loR),
                    datalogger.createCV("first", first),
                    datalogger.createCV("last", last),
                    datalogger.createCV("semiPeriod", this.semiPeriod))
            }
        }


        // Transform a point on the off-centre projected Ellipse back onto the centred Spin-Circle 
        // and return its angle (in radians anticlockwise from the horizontal U-axis)
        project(uRaw: number, vRaw: number): number {
            // shift the point to give a vector from the origin
            let u = uRaw - this.uOff
            let v = vRaw - this.vOff
            // rotate clockwise, using the INVERSE of the major-axis angle, theta
            let uNew = u * this.cosTheta - v * this.sinTheta
            let vNew = u * this.sinTheta + v * this.cosTheta
            // scale up V to match U, balancing the axes and making the Ellipse circular
            let vScaled = vNew * this.scale
            // derive projected angle (undoing the rotation we just applied)
            let angle = Math.atan2(vScaled, uNew) + this.theta
            if (logging) {
                datalogger.setColumnTitles("uRaw", "vRaw", "u", "v", "uNew", "vNew", "vScaled", "angle")
                datalogger.log(datalogger.createCV("uRaw", uRaw),
                    datalogger.createCV("vRaw", vRaw),
                    datalogger.createCV("u", u),
                    datalogger.createCV("v", v),
                    datalogger.createCV("uNew", uNew),
                    datalogger.createCV("vNew", vNew),
                    datalogger.createCV("vScaled", vScaled),
                    datalogger.createCV("angle", angle))
            }

            return angle
        }
        
    }

    // GLOBALS
   
    let scanTimes: number[]  = [] // sequence of time-stamps for scanned readings 
    let scanData: number[][] = [] // scanned sequence of [X,Y,Z] magnetometer readings
    let views: Ellipse[] = [] // the three possible elliptical views of the Spin-Circle
    let bestView = -1
    let uDim = -1 // the best "horizontal" axis (pointing East) for transformed readings (called U)
    let vDim = -1 // the best "vertical" axis (pointing North) for transformed readings (called V)
    let toNorth = 0 // the angular bias to be added (so that North = 0)
    let strength = 0 // the average magnetic field-strength observed by the magnetometer
    let period = 0 // average rotation time derived from scanData[]
    //let uFlip = 1 // set to -1 if uDim polarity is inverted
    //let vFlip = 1 // set to -1 if vDim polarity is inverted

    let logging = true  // logging mode flag
    let testing = false  // test mode flag
    let test = 0 // selector for test cases
    let uData: number[] = [] // U values for test cases
    let vData: number[] = [] // V values for test cases

    // EXPORTED USER INTERFACES   

    /** 
     * Assuming the buggy is currently spinning clockwise on the spot, capture a 
     * time-stamped sequence of magnetometer readings from which to set up the compass.
     * 
     * Whwn complete, analyse the scanned data to prepare for reading compass-headings.
     * Returns zero if successful, or a negative error code:
     *
     *      -1 : NOT ENOUGH SCAN DATA
     *
     *      -2 : FIELD STRENGTH TOO WEAK
     *
     *      -3 : NOT ENOUGH SCAN ROTATION
     * 
     * @param ms scanning-time in millisecs (long enough for more than one full rotation)    
     */

    //% block="scan for (ms) $ms" 
    //% inlineInputMode=inline 
    //% ms.shadow="timePicker" 
    //% ms.defl=0 
    //% weight=90 

    export function scan(ms: number): number {
        // Every ~30 ms over the specified duration (generally a couple of seconds),
        // magnetometer readings are sampled and a new [X,Y,Z] triple added to the scanData[] array.
        // A timestamp for each sample is also recorded in the scanTimes[] array.
        if (logging) {
            datalogger.deleteLog()
            datalogger.includeTimestamp(FlashLogTimeStampFormat.None)
            //datalogger.setColumnTitles("trace")
            //datalogger.log(datalogger.createCV("trace", 1))
        }
        if (testing) {
            simulateScan("Y-up") // Y-Axis vertical; spun around Y-Axis
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

    // Now analyse the scan-data to decide how best to use the magnetometer readings.

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
        // To assess the range and offset, first find the raw extremes.
        for (let i = 0; i < nSamples; i++) {
            xhi = Math.max(xhi, scanData[i][Dim.X])
            yhi = Math.max(yhi, scanData[i][Dim.Y])
            zhi = Math.max(zhi, scanData[i][Dim.Z])
            xlo = Math.min(xlo, scanData[i][Dim.X])
            ylo = Math.min(ylo, scanData[i][Dim.Y])
            zlo = Math.min(zlo, scanData[i][Dim.Z])
        }

        // Use the mean of these extremes as normalisation offsets
        let xOff = (xhi + xlo) / 2
        let yOff = (yhi + ylo) / 2
        let zOff = (zhi + zlo) / 2

        // create an Ellipse instance for analysing each possible view
        views.push(new Ellipse("XY", Dim.X, Dim.Y, xOff, yOff))
        views.push(new Ellipse("YZ", Dim.Y, Dim.Z, yOff, zOff))
        views.push(new Ellipse("ZX", Dim.Z, Dim.X, zOff, xOff))

        // Now assess eccentricities by looking for the shortest and longest radii for each Ellipse.
        strength = 0 // initialise global field-strength-squared accumulator
        
        // initialise the smoothed radii and their limits to their first (unsmoothed) values
        let xRaw = scanData[0][Dim.X]
        let yRaw = scanData[0][Dim.Y]
        let zRaw = scanData[0][Dim.Z]

        views[View.XY].firstSample(xRaw, yRaw)
        views[View.YZ].firstSample(yRaw, zRaw)
        views[View.ZX].firstSample(zRaw, xRaw)

        if (logging) {
            // prepare for logging peaks within Ellipse.nextSample() function...
            datalogger.setColumnTitles("view", "index", "loRsq", "hiRsq")
        }

        // process scan from 2nd sample onwards...
        for (let j = 1; j < nSamples; j++) {
            xRaw = scanData[j][Dim.X]
            yRaw = scanData[j][Dim.Y]
            zRaw = scanData[j][Dim.Z]
            views[View.XY].nextSample(j, xRaw, yRaw) // projection in XY plane
            views[View.YZ].nextSample(j, yRaw, zRaw) // projection in YZ plane
            views[View.ZX].nextSample(j, zRaw, xRaw) // projection in ZX plane

            // Square of overall 3-D field strength is sum of squares in each View
            // Here we will end up accumulating each one twice over!
            strength += views[View.XY].rSq
            strength += views[View.YZ].rSq
            strength += views[View.ZX].rSq
        }

   
        // check average overall field-strength (undoing the double-counting)
        strength = Math.sqrt((strength / 2) / nSamples)
        if (strength < MarginalField) {
            return -2  // "FIELD STRENGTH TOO WEAK"
        }

        // ?? check that we have collected enough peaks in each View...
        if ((views[View.XY].peaks.length < 3) 
          ||(views[View.YZ].peaks.length < 3)
          ||(views[View.ZX].peaks.length < 3) ) {
            return -3 // "NOT ENOUGH SCAN ROTATION"
        }


        // process the detected Ellipse axes
        views[View.XY].analyse()
        views[View.YZ].analyse()
        views[View.ZX].analyse()


        // Choose the "roundest" Ellipse, the one with lowest (eccentricity = scale).
        // While we're at it, derive the rotation period, in case anyone's interested.
        // The other two (more eccentric) Ellipses will be more reliable for deriving period,
        // so we collect all three semiPeriods, but then subtract the [bestView] version.

        bestView = View.XY
        period = -views[View.XY].semiPeriod
        
        if (views[View.YZ].scale < views[bestView].scale) { 
            bestView = View.YZ
            period = -views[View.YZ].semiPeriod
        }

        if (views[View.ZX].scale < views[bestView].scale) {
            bestView = View.YZ
            period = -views[View.ZX].semiPeriod
        }

        period += views[View.XY].semiPeriod 
                + views[View.YZ].semiPeriod 
                + views[View.ZX].semiPeriod

        basic.clearScreen()
        basic.pause(100)
        basic.showString(views[bestView].plane)
        basic.pause(300)

        // extract into globals for future brevity and efficiency...
        uDim = views[bestView].uDim
        vDim = views[bestView].vDim

        // we've now finished with scanning data, so release its memory
        scanTimes = []
        scanData = []

        // SUCCESS!
        return 0
    }



    /**
     * Read the magnetometer and register the buggy's current direction as "North",
     * (i.e. the direction that will in future return zero as its heading).
     * The actual direction the buggy is pointing when this function is called could be
     * Magnetic North; True North (compensating for declination); or any convenient
     * direction from which to measure subsequent heading angles.
     */
    //% block="set North" 
    //% inlineInputMode=inline 
    //% weight=80 
    export function setNorth() {
        // We have successfully set up the projection parameters. Now we need to relate them to North.
        // Take the average of seven new readings to get a stable fix on the current heading
        let uRaw = 0
        let vRaw = 0
        if (testing) { //arbitrarily choose 10th test reading for North
            uRaw = uData[0]
            vRaw = vData[0]
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

        

        // get the projection angle of the [uRaw,vRaw] vector on the Spin-Circle 
        // and remember this as the (global) fixed bias to North
        toNorth = views[bestView].project(uRaw, vRaw)

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
                datalogger.createCV("uOff", views[bestView].uOff),
                datalogger.createCV("vOff", views[bestView].vOff),
                datalogger.createCV("theta", views[bestView].theta),
                datalogger.createCV("scale", views[bestView].scale),
                datalogger.createCV("period", period),
                datalogger.createCV("toNorth", toNorth),
                datalogger.createCV("strength", strength))
        }
    }
  



    /**
     * Read the magnetometer and return the current heading of the buggy in degrees
     */
    //% block="degrees" 
    //% inlineInputMode=inline 
    //% weight=70
    export function degrees(): number {
        // read the magnetometer (seven times) and return the current heading in degrees from North
        let uRaw = 0
        let vRaw = 0
        if (testing) {
            uRaw = uData[test]
            vRaw = vData[test]
            test = (test+1) % uData.length
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

        // project reading from Ellipse-view to Spin-Circle, as radians anticlockwise from U-axis
        let onCircle = views[bestView].project(uRaw, vRaw)

        // subtracting "toNorth" gives an angle measured anticlockwise from North;
        // negating this gives an angle measured clockwise from North
        let angle = radians2degrees(toNorth - onCircle)
        return angle
    }

    /**
     * Return average RPM of the most recent scan 
     */
    //% block="scanRPM" 
    //% inlineInputMode=inline 
    //% weight=50 
    export function scanRPM(): number {
        if (period <= 0) {
            return -4 // ERROR: SUCCESSFUL SCAN IS NEEDED FIRST
        } else {
            return 60000 / period
        }
    }


    //% block="set test mode: $turnOn"
    //% inlineInputMode=inline 
    //% weight=40
    export function setTestMode(turnOn: boolean) {
        testing = turnOn
    }

    //% block="set logging mode: $loggingOn"
    //% inlineInputMode=inline 
    //% weight=30
    export function setLogMode(loggingOn: boolean) {
        logging = loggingOn
    }

    export function isLogging(): boolean {
        return logging
    }



// UTILITY FUNCTIONS

    // convert angle from radians to degrees in the range 0..360
    function radians2degrees(angle: number): number {
        let degrees = ((57.29578 * angle) + 360) % 360
        return degrees
    }

    // While debugging, it is useful to re-use predictable sample data for a variety of use-cases 
    // that has been captured from live runs using datalogger.
    function simulateScan(dataset: string) {
        let xData: number[] = []
        let yData: number[] = []
        let zData: number[] = []
        switch (dataset) {
            case "Y-up": // Y-axis vertical; spun around Y-axis
                scanTimes = [69597, 69625, 69653, 69681, 69709, 69737, 69765, 69793, 69821, 69851, 69877, 69905, 69933, 69961, 69989, 70017, 70045, 70073, 70101, 70129, 70157, 70185, 70213, 70241, 70269, 70297, 70325, 70353, 70381, 70412, 70437, 70465, 70493, 70521, 70549, 70577, 70605, 70633, 70661, 70689, 70717, 70745, 70773, 70801, 70829, 70857, 70888, 70913, 70941, 70969, 70997, 71025, 71053, 71081, 71109, 71138, 71165, 71193, 71221, 71250, 71278, 71305, 71336, 71361, 71389, 71417, 71445, 71473, 71501, 71529, 71557, 71585, 71614, 71642, 71670, 71701, 71729, 71757, 71785, 71813, 71841, 71870, 71898, 71926, 71954, 71982, 72010, 72038, 72069, 72097, 72126, 72153, 72181, 72209, 72238, 72266, 72294, 72322, 72350, 72378, 72409, 72437, 72465, 72494, 72521, 72549, 72578, 72609, 72637, 72665, 72693, 72721, 72750, 72778, 72806, 72834, 72862, 72893, 72921, 72949, 72978, 73006, 73034, 73062, 73090, 73121, 73149, 73177, 73206, 73234, 73262, 73290, 73321, 73349, 73377, 73406, 73434, 73462, 73490, 73521, 73550, 73578, 73606, 73634, 73662, 73693, 73721, 73750, 73778, 73806, 73837, 73865, 73894, 73922, 73950, 73981, 74010, 74038, 74066, 74097, 74125, 74154, 74182, 74213, 74241, 74270, 74298, 74329, 74358, 74386, 74414, 74445, 74474, 74502, 74533, 74562, 74593, 74626, 74653, 74681, 74709, 74743, 74770, 74798, 74826, 74854, 74882, 74910, 74938, 74966, 74994, 75022, 75050, 75078, 75106, 75134, 75162, 75190, 75218, 75246, 75274, 75302, 75330, 75358, 75386, 75420, 75446, 75474, 75502, 75530, 75558, 75586, 75614, 75642, 75670, 75698, 75726, 75754, 75782, 75810, 75838, 75866, 75894, 75922, 75950, 75978, 76012, 76038, 76066, 76094, 76122, 76150, 76178, 76206, 76234, 76262, 76290, 76318, 76346, 76374, 76402, 76430, 76458, 76486, 76514, 76549, 76578, 76606, 76634, 76662, 76690, 76718, 76746, 76774, 76802, 76830, 76858, 76886, 76914, 76942, 76971, 77005, 77034, 77062, 77090, 77118, 77146, 77174, 77202, 77230, 77258, 77286, 77314, 77342, 77370, 77398, 77433, 77462, 77490, 77518, 77546, 77574, 77602, 77630, 77658, 77686, 77714, 77742, 77770, 77798, 77833, 77862, 77890, 77918, 77953, 77982, 78010, 78038, 78066, 78094, 78122, 78150, 78178, 78206, 78234, 78269, 78298, 78326, 78354, 78382, 78410, 78438, 78466, 78494, 78522, 78557, 78586, 78614, 78642, 78670, 78698, 78726, 78754, 78789, 78818, 78846, 78874, 78902, 78930, 78958, 78993, 79022, 79050, 79078, 79106, 79134, 79169, 79198, 79226, 79254, 79282, 79310, 79345, 79374, 79402, 79430, 79458, 79493, 79522, 79550, 79578, 79606, 79642, 79670, 79698, 79726, 79754, 79789, 79818, 79846, 79874, 79909, 79938, 79966, 79995, 80029, 80058, 80086, 80122, 80150, 80178, 80214, 80242, 80270, 80305, 80333, 80362, 80389, 80417, 80445, 80473, 80501, 80529, 80557, 80585, 80613, 80651, 80678, 80706, 80734, 80762, 80790, 80818, 80846, 80874, 80902, 80930, 80958, 80986, 81014, 81042, 81070, 81098, 81126, 81154, 81182, 81210, 81238, 81266, 81294, 81333, 81362, 81390, 81418, 81446, 81474, 81502, 81530, 81558, 81586, 81614, 81642, 81670, 81698, 81726, 81754, 81782, 81810, 81838, 81866, 81895, 81922, 81961, 81990, 82018, 82046, 82074, 82102, 82130, 82158, 82186, 82214, 82242, 82271, 82299, 82326, 82355, 82383, 82411, 82439, 82467, 82506, 82534, 82562, 82590, 82618, 82646, 82674, 82702, 82730, 82759, 82787, 82815, 82843, 82871, 82899, 82927, 82955, 82994, 83022, 83050, 83078, 83106, 83135, 83162, 83190, 83218, 83247, 83275, 83303, 83331, 83359, 83387, 83415, 83454, 83482, 83510, 83538, 83566, 83606, 83633, 83661, 83702, 83730, 83758, 83786, 83814, 83842, 83871, 83898, 83926, 83954, 83982, 84010, 84038, 84066, 84094, 84123, 84150, 84179, 84207, 84235, 84263, 84291, 84319, 84347, 84375, 84403, 84431]
                xData = [-46.35, -46.5, -46.8, -47.4, -47.7, -48, -48.45, -49.05, -49.2, -49.95, -50.4, -51, -52.35, -53.85, -55.35, -57.15, -58.2, -59.85, -62.1, -63.15, -64.35, -65.4, -66.75, -69.15, -70.65, -71.85, -73.05, -75.15, -77.4, -79.65, -81.3, -82.8, -84.6, -87.6, -89.55, -90.6, -93, -95.25, -98.1, -101.25, -103.65, -106.8, -111.15, -114.3, -118.05, -121.35, -124.2, -126.45, -130.5, -132.9, -135.3, -136.95, -139.2, -141.45, -145.05, -146.1, -148.5, -150.6, -153.6, -157.65, -161.85, -164.85, -168.45, -172.5, -176.85, -181.05, -183.9, -186.75, -190.8, -195.6, -200.1, -204.15, -208.05, -212.1, -216, -219.9, -224.1, -227.4, -230.7, -234.45, -238.05, -241.65, -245.4, -248.7, -252.6, -256.65, -259.95, -263.25, -266.7, -269.25, -270.3, -271.5, -273.3, -275.4, -277.8, -279.45, -281.55, -284.55, -286.2, -287.7, -289.2, -290.4, -290.85, -291.6, -292.05, -293.25, -293.7, -294, -295.05, -296.4, -297, -297.9, -298.35, -299.4, -300.6, -301.2, -301.35, -301.5, -302.4, -302.85, -303.15, -302.55, -302.55, -304.65, -305.55, -305.85, -306, -306.75, -307.35, -307.2, -306.15, -305.1, -304.35, -303.9, -301.95, -301.2, -300.15, -298.05, -297, -295.8, -294.6, -294.3, -292.8, -291.45, -290.7, -289.65, -288, -286.5, -285, -284.1, -281.55, -278.85, -276, -273.3, -270.45, -267.6, -263.85, -261.45, -259.5, -257.1, -254.4, -252.3, -249.9, -247.5, -245.55, -243.3, -242.1, -240.9, -239.4, -237.6, -235.35, -233.25, -230.55, -227.1, -223.95, -220.95, -217.35, -215.4, -212.25, -209.7, -206.85, -205.05, -202.35, -200.1, -197.25, -194.85, -192.45, -190.05, -187.5, -184.8, -181.95, -178.8, -177.3, -175.2, -173.25, -169.65, -167.25, -164.55, -162.15, -157.8, -154.95, -151.65, -149.7, -146.1, -143.4, -141.15, -139.5, -136.05, -134.25, -131.25, -129.3, -126.6, -123.75, -120.3, -118.65, -115.5, -112.8, -111.45, -109.8, -108.3, -107.25, -104.85, -103.35, -100.95, -98.25, -96.45, -93.9, -92.1, -90.15, -88.8, -88.05, -86.25, -85.05, -84.15, -82.8, -82.5, -80.1, -78.75, -77.7, -75.6, -74.4, -72.6, -70.5, -70.05, -69, -66.75, -65.55, -64.2, -63.3, -61.8, -60.3, -59.55, -59.4, -58.2, -56.85, -56.4, -56.1, -55.5, -53.7, -53.55, -53.4, -52.95, -51.6, -50.25, -49.95, -49.35, -48.45, -48.15, -48.15, -48.15, -48.45, -47.55, -46.95, -46.05, -44.85, -43.2, -41.7, -41.85, -41.4, -42, -42.45, -42.75, -43.05, -44.4, -43.35, -44.1, -44.1, -43.8, -44.4, -45.45, -46.2, -46.65, -46.5, -46.95, -48.45, -49.05, -50.1, -49.2, -50.55, -51.6, -52.5, -52.05, -52.2, -51.6, -52.35, -52.35, -52.05, -52.8, -52.95, -53.55, -54.45, -55.05, -55.5, -57.15, -56.85, -57.6, -57.75, -58.2, -58.95, -59.1, -58.65, -59.85, -61.05, -61.95, -63, -63.45, -65.25, -66.3, -66, -66.45, -67.2, -67.05, -67.65, -67.5, -68.4, -69.45, -70.35, -70.5, -71.55, -72.75, -73.95, -74.1, -74.7, -75.45, -76.95, -78.75, -79.8, -80.55, -82.5, -84.15, -85.95, -88.05, -88.8, -89.7, -91.65, -93.15, -94.5, -95.85, -97.95, -98.55, -101.55, -103.2, -104.1, -106.5, -108.3, -108.75, -111, -111.9, -113.55, -115.35, -117, -119.25, -121.5, -124.8, -127.65, -128.85, -131.1, -132.75, -134.55, -136.95, -139.2, -140.55, -143.55, -146.85, -149.55, -151.5, -153.3, -155.7, -157.5, -160.05, -162.15, -165.15, -168.3, -171.45, -174.3, -177.6, -180.9, -183.6, -185.85, -187.95, -190.2, -191.55, -194.7, -197.1, -200.1, -202.35, -204.45, -207.45, -209.4, -210.45, -211.5, -213.9, -216, -218.25, -220.2, -222.45, -226.2, -229.35, -231.75, -234, -237.45, -239.85, -243.15, -244.5, -246.9, -248.25, -251.1, -252.75, -254.55, -256.5, -258.15, -260.55, -262.35, -264.75, -266.85, -268.5, -269.25, -271.95, -272.4, -273.6, -274.65, -274.35, -275.1, -277.05, -277.8, -279.75, -281.1, -280.8, -282.45, -283.65, -284.7, -285.3, -285.15, -285.75, -286.65, -287.55, -289.35, -288.75, -288.9, -289.65, -291, -292.35, -292.65, -291.75, -292.8, -294.15, -295.35, -295.05, -295.65, -296.55, -298.65, -299.4, -300.15, -300.3, -301.05, -301.05, -301.5, -299.7, -300, -300.6, -300.3, -300, -300.6, -301.5, -302.7, -303.9, -303, -304.2, -305.25, -305.1, -304.95, -305.25, -304.8, -304.95, -304.65, -304.95, -304.95, -304.8, -304.35, -304.35, -304.5, -304.05, -303.3, -303.9, -303.75, -303.45, -303.6, -303.15, -303.45, -303.6, -301.8, -301.8, -301.5, -301.5, -301.65, -300.9, -299.55, -300.3, -298.5]
                yData = [-33.3, -32.85, -32.55, -32.1, -31.95, -31.95, -33, -33, -33.9, -33.6, -34.05, -34.5, -34.8, -34.95, -33.75, -32.7, -32.85, -32.7, -32.55, -32.1, -31.35, -32.4, -32.7, -33.45, -33.9, -34.35, -34.65, -35.1, -35.1, -35.4, -34.95, -34.8, -33.9, -34.2, -34.2, -34.05, -34.05, -34.05, -34.05, -34.35, -33.6, -33.3, -32.7, -32.85, -32.7, -33.45, -34.05, -35.1, -35.4, -34.95, -34.5, -34.8, -33, -32.1, -31.2, -30.75, -31.65, -31.65, -31.05, -32.4, -32.85, -32.25, -32.7, -32.7, -32.4, -32.4, -31.95, -31.8, -32.7, -32.55, -32.1, -31.95, -31.8, -31.8, -31.35, -30.9, -30.3, -30.45, -30.9, -31.2, -30.45, -30.9, -30.6, -31.2, -31.5, -31.2, -30.9, -30.45, -30, -30, -29.7, -28.95, -28.35, -28.65, -28.65, -28.35, -28.5, -28.2, -28.35, -28.95, -29.7, -30.3, -30.9, -31.2, -31.8, -32.7, -32.25, -31.5, -31.35, -31.2, -31.2, -31.2, -31.2, -31.65, -31.2, -31.35, -30.75, -30.45, -29.7, -29.1, -29.25, -29.7, -29.85, -29.85, -30.3, -30.15, -30, -29.1, -28.8, -28.5, -29.1, -27.9, -28.2, -28.2, -28.5, -28.95, -28.5, -27.75, -28.35, -28.65, -28.5, -28.8, -28.5, -28.65, -29.25, -28.65, -28.05, -28.2, -28.65, -27.6, -28.05, -27.15, -27.45, -27.75, -26.7, -25.2, -25.35, -24.3, -25.95, -25.5, -24.6, -25.5, -25.2, -25.35, -26.55, -25.65, -26.1, -27.3, -27.3, -27.6, -27.45, -27.15, -26.55, -26.7, -26.4, -26.1, -26.55, -27.3, -27, -26.85, -26.55, -26.7, -27, -27.3, -26.85, -27.3, -28.5, -29.1, -29.25, -28.65, -29.4, -29.55, -29.85, -28.65, -29.1, -28.95, -29.85, -30, -29.85, -30.15, -30.45, -29.55, -28.95, -29.1, -28.2, -29.25, -27.6, -27.3, -27.45, -28.5, -29.1, -29.1, -28.5, -30, -30.45, -31.2, -31.05, -30.15, -30.45, -30.75, -31.05, -31.65, -30.75, -30.9, -31.2, -30.6, -30.9, -30.3, -30, -30.6, -30, -30.15, -30.75, -30.15, -29.7, -29.85, -30.6, -30.3, -30.9, -30.3, -30.75, -31.35, -31.8, -31.2, -31.2, -30.9, -32.25, -31.95, -31.8, -31.95, -32.7, -33.15, -32.4, -32.85, -32.55, -33.45, -33.45, -32.85, -33.15, -33.15, -32.4, -32.85, -32.4, -31.35, -31.2, -32.25, -32.25, -33.15, -33.9, -33.45, -34.35, -34.8, -33.6, -34.35, -33.15, -32.4, -32.7, -32.7, -31.8, -32.25, -31.35, -32.1, -32.85, -32.4, -32.55, -33.6, -33.9, -34.95, -34.95, -34.65, -34.8, -34.65, -34.35, -34.05, -33.6, -33.6, -34.05, -34.05, -33.75, -34.05, -33.6, -34.2, -33.75, -33.15, -33.75, -33.75, -33.45, -33.45, -33.15, -33.75, -34.2, -34.05, -34.65, -34.5, -34.95, -34.8, -34.65, -35.1, -34.95, -34.65, -34.95, -34.8, -35.1, -34.2, -33.15, -33.75, -33.75, -34.05, -33.9, -33.75, -34.2, -34.65, -33.75, -33.75, -32.55, -33.45, -32.25, -32.4, -31.8, -32.4, -31.65, -31.35, -30.9, -31.35, -30.9, -30.75, -30.45, -30.15, -30.45, -30, -30.9, -31.35, -31.65, -31.5, -33.15, -34.05, -33.9, -32.85, -33, -33, -33, -31.95, -31.8, -32.4, -32.7, -31.8, -31.05, -31.35, -31.05, -31.2, -31.2, -31.65, -32.7, -33, -33.6, -34.35, -34.2, -33.45, -33.6, -32.7, -32.25, -31.35, -31.05, -30.45, -30.9, -30.15, -30.3, -31.35, -30.9, -31.35, -31.35, -31.35, -32.4, -32.7, -31.95, -31.95, -31.2, -31.2, -30.6, -29.55, -28.65, -28.65, -29.55, -29.7, -30.15, -30.3, -30.75, -31.2, -32.25, -31.8, -31.65, -30.9, -31.35, -31.35, -31.65, -30.6, -30.6, -31.05, -31.8, -31.65, -31.05, -30.45, -30.3, -30.3, -29.55, -29.25, -29.25, -29.25, -30, -30.15, -29.85, -30.45, -30.6, -30, -30.6, -30.75, -30.6, -31.05, -31.35, -31.8, -32.25, -32.25, -31.8, -32.25, -31.5, -30.9, -29.4, -29.25, -29.1, -29.4, -29.7, -30.6, -30.9, -31.35, -31.5, -31.05, -30.9, -30.6, -29.85, -29.7, -29.85, -30.15, -30.15, -29.55, -29.1, -29.25, -29.1, -29.7, -29.1, -29.25, -29.55, -29.55, -29.55, -29.7, -28.8, -28.2, -28.65, -29.1, -29.85, -30.3, -30.45, -30.3, -30.9, -29.55, -29.1, -28.95, -28.5, -28.2, -27.9, -27.75, -28.8, -28.5, -27.75, -27.15, -27.45, -27.3, -27.6, -27.15, -27.6, -27.75, -28.2, -27.6, -28.05, -28.05, -27.75, -27, -26.85, -27.45, -28.35, -28.8, -28.05, -28.65, -29.1, -29.55, -28.5, -27.9, -27.9, -27.45, -26.2]
                zData = [-71.1, -74.25, -77.1, -80.1, -84, -86.85, -88.5, -92.7, -95.1, -97.35, -99.3, -102.3, -106.65, -109.95, -110.55, -112.8, -115.8, -119.25, -121.35, -123, -125.7, -128.85, -131.4, -133.65, -135.6, -138.15, -140.1, -142.35, -144.3, -146.4, -148.65, -149.7, -150.6, -153.45, -155.55, -157.05, -159.15, -159.45, -162.75, -165.6, -166.95, -168.9, -171.45, -173.25, -176.7, -177.75, -178.5, -179.4, -180.15, -180.75, -182.25, -182.85, -183.9, -184.65, -184.95, -185.25, -187.5, -187.95, -187.95, -187.95, -189, -190.35, -190.8, -188.7, -188.4, -188.25, -188.4, -188.1, -186.9, -186.75, -186.75, -184.8, -184.8, -184.35, -182.85, -182.25, -180.3, -179.25, -179.55, -177.45, -174.15, -172.35, -169.95, -168, -165.45, -162.3, -159.75, -158.25, -155.85, -153.15, -150.9, -149.4, -146.4, -144.3, -141.75, -139.2, -137.4, -135, -131.85, -130.35, -128.7, -127.05, -124.8, -121.95, -121.2, -118.95, -116.7, -113.4, -111.15, -109.35, -107.85, -103.65, -101.55, -98.7, -96.6, -93.15, -90.6, -88.2, -86.55, -83.7, -80.55, -77.85, -75.3, -71.25, -67.8, -63.9, -61.65, -59.4, -55.8, -52.5, -51, -47.7, -44.55, -40.95, -37.65, -35.1, -32.1, -27.15, -23.25, -20.4, -15.75, -11.7, -8.55, -5.85, -4.05, -2.25, 1.05, 2.85, 4.8, 8.25, 11.4, 14.7, 16.95, 19.35, 21.3, 24.15, 25.65, 27.9, 30.15, 33.3, 34.5, 36.6, 38.4, 40.2, 41.4, 41.85, 41.85, 43.35, 45.45, 45.15, 45.6, 46.65, 47.85, 48.6, 48.75, 48.75, 50.85, 52.65, 53.25, 54.3, 55.95, 57.75, 58.8, 59.1, 59.4, 60.3, 59.7, 59.25, 60.75, 60.45, 59.55, 60.3, 60.15, 61.35, 61.8, 60.3, 61.2, 61.65, 60.75, 60.9, 60, 60.15, 60, 58.5, 58.65, 58.5, 58.2, 57.75, 56.55, 54.9, 54.3, 52.95, 52.05, 49.95, 49.8, 49.2, 48.75, 46.8, 46.05, 44.1, 43.2, 40.95, 39.3, 39.15, 38.7, 36.9, 36.15, 35.55, 34.05, 33.75, 30.9, 29.55, 28.5, 27.6, 24.75, 24.15, 22.5, 21, 19.95, 18.3, 16.5, 15.9, 13.65, 11.55, 9.75, 8.7, 6.6, 4.35, 1.65, 0.9, -0.6, -1.95, -5.85, -6.9, -7.8, -9.45, -12.15, -14.4, -15.75, -16.8, -18.75, -21.75, -22.95, -24.75, -26.7, -29.7, -31.05, -33.3, -34.5, -35.25, -37.2, -38.4, -39.15, -40.8, -43.35, -45.45, -48.15, -50.4, -53.85, -56.1, -58.05, -58.95, -61.95, -64.5, -66.15, -66.75, -69.6, -73.05, -75.6, -76.8, -78.9, -81.75, -84, -85.95, -86.85, -88.8, -91.05, -92.55, -93.45, -94.65, -95.7, -97.35, -97.65, -98.7, -100.95, -102.3, -103.65, -104.4, -105.75, -107.85, -108.15, -107.4, -108.75, -110.4, -111.45, -112.2, -112.65, -114.75, -116.7, -116.7, -116.85, -118.35, -120, -121.8, -122.85, -125.55, -128.25, -129.15, -129.45, -131.4, -132.45, -132.9, -132, -131.85, -135.3, -136.8, -136.65, -137.55, -139.5, -141.45, -142.35, -142.2, -143.4, -144.9, -145.65, -145.5, -145.95, -147.9, -147.9, -148.95, -150, -152.25, -153.9, -155.7, -157.05, -160.5, -162.3, -162.6, -163.5, -165.75, -166.65, -168, -168.15, -170.25, -171.9, -172.5, -172.5, -174, -174.75, -174.9, -175.2, -175.95, -177.3, -178.8, -179.1, -180.15, -181.5, -181.5, -181.5, -182.1, -183.15, -184.05, -183.6, -184.2, -185.4, -187.2, -188.1, -187.5, -189, -189.9, -190.95, -191.1, -190.65, -190.65, -191.7, -190.65, -190.5, -190.8, -190.2, -189.9, -189.3, -188.55, -187.95, -188.4, -186.15, -185.4, -185.25, -184.65, -183.3, -182.55, -181.95, -181.95, -181.2, -180.15, -179.7, -179.4, -177.6, -175.65, -175.8, -175.65, -174.6, -173.55, -172.05, -171.9, -172.35, -169.35, -167.85, -166.95, -165.15, -163.95, -162.45, -159.75, -160.2, -158.4, -156.75, -155.1, -154.2, -152.55, -150.3, -146.85, -146.1, -145.35, -142.95, -140.7, -140.25, -140.1, -139.65, -137.4, -135.3, -136.35, -135.75, -134.1, -132.3, -131.25, -130.2, -129.3, -126.45, -125.7, -124.2, -122.1, -120.3, -118.95, -118.05, -117.6, -115.35, -113.85, -113.25, -112.65, -110.85, -109.05, -107.1, -105.9, -103.65, -101.85, -99.15, -97.95, -96.3, -94.2, -92.25, -91.2, -89.4, -88.05, -87.45, -85.95, -85.5, -83.7, -82.8, -82.05, -82.2, -79.5, -77.85, -76.2, -75, -72.3, -70.65, -67.95, -66.9, -63.45, -60.6, -58.5, -57.45, -55.05, -52.35, -49.65, -49.65, -48.15, -45.75, -42.6, -40.5, -38.25, -36.15, -33.45, -31.05, -30.3, -28.95, -26.4, -24.75, -21.9, -20.25]
                uData = [-34.05, -32.85, -31.2, -24.6, -29.1, -32.4, -34.2, -33, -34.05]
                vData = [-99.3, -189, -107.85, 34.5, 54.3, -18.75, -108.15, -166.65, -99.3]
                break
            default: // tetrahedral: spin-axis 45 degrees to X, Y & Z
                scanTimes = [51901, 51929, 51957, 51985, 52013, 52041, 52069, 52097, 52125, 52155, 52182, 52209, 52237, 52265, 52293, 52321, 52349, 52377, 52405, 52433, 52461, 52489, 52517, 52545, 52573, 52601, 52629, 52657, 52685, 52716, 52741, 52769, 52797, 52825, 52853, 52881, 52909, 52937, 52965, 52993, 53021, 53049, 53077, 53105, 53134, 53161, 53192, 53217, 53245, 53273, 53301, 53329, 53357, 53385, 53413, 53441, 53469, 53497, 53525, 53553, 53581, 53609, 53640, 53669, 53697, 53725, 53754, 53781, 53809, 53837, 53865, 53893, 53921, 53950, 53978, 54009, 54037, 54065, 54093, 54121, 54149, 54177, 54205, 54233, 54262, 54289, 54318, 54346, 54377, 54405, 54433, 54461, 54489, 54517, 54545, 54573, 54601, 54630, 54658, 54686, 54716, 54741, 54769, 54797, 54825, 54853, 54882, 54913, 54941, 54969, 54997, 55025, 55054, 55081, 55110, 55138, 55166, 55196, 55221, 55249, 55277, 55306, 55334, 55362, 55390, 55421, 55449, 55477, 55505, 55533, 55562, 55590, 55621, 55649, 55677, 55705, 55734, 55762, 55790, 55821, 55849, 55877, 55906, 55934, 55962, 55993, 56021, 56049, 56077, 56106, 56137, 56165, 56193, 56221, 56250, 56278, 56309, 56337, 56366, 56394, 56425, 56453, 56482, 56510, 56541, 56569, 56598, 56626, 56657, 56686, 56714, 56742, 56773, 56802, 56830, 56861, 56890, 56921, 56949, 56977, 57005, 57038, 57065, 57094, 57122, 57150, 57178, 57206, 57234, 57262, 57290, 57318, 57346, 57374, 57402, 57430, 57458, 57486, 57514, 57542, 57570, 57598, 57626, 57654, 57682, 57710, 57744, 57770, 57798, 57826, 57854, 57882, 57910, 57938, 57966, 57994, 58022, 58050, 58078, 58106, 58134, 58162, 58190, 58218, 58246, 58274, 58308, 58333, 58362, 58390, 58418, 58446, 58474, 58502, 58530, 58558, 58586, 58614, 58642, 58670, 58698, 58726, 58754, 58782, 58810, 58845, 58874, 58902, 58930, 58958, 58986, 59014, 59042, 59070, 59098, 59126, 59154, 59182, 59210, 59238, 59266, 59294, 59329, 59358, 59386, 59414, 59442, 59470, 59498, 59526, 59554, 59582, 59610, 59638, 59666, 59694, 59722, 59757, 59786, 59814, 59842, 59870, 59898, 59926, 59954, 59982, 60010, 60038, 60066, 60094, 60122, 60150, 60185, 60214, 60249, 60278, 60306, 60334, 60362, 60390, 60418, 60446, 60474, 60502, 60530, 60558, 60593, 60621, 60650, 60678, 60706, 60734, 60762, 60790, 60818, 60853, 60881, 60910, 60938, 60966, 60994, 61022, 61050, 61078, 61113, 61142, 61170, 61198, 61226, 61254, 61282, 61310, 61345, 61374, 61402, 61430, 61458, 61486, 61514, 61549, 61578, 61606, 61634, 61662, 61690, 61725, 61754, 61782, 61810, 61838, 61873, 61902, 61930, 61958, 61986, 62021, 62050, 62078, 62106, 62134, 62169, 62198, 62226, 62254, 62289, 62318, 62346, 62375, 62409, 62438, 62466, 62502, 62530, 62558, 62594, 62622, 62650, 62685, 62713, 62741, 62769, 62797, 62825, 62853, 62881, 62909, 62938, 62965, 62993, 63031, 63058, 63086, 63114, 63142, 63171, 63198, 63226, 63254, 63282, 63310, 63338, 63366, 63394, 63422, 63450, 63478, 63506, 63534, 63562, 63590, 63618, 63646, 63674, 63712, 63738, 63766, 63794, 63822, 63851, 63878, 63906, 63934, 63962, 63990, 64018, 64046, 64074, 64103, 64130, 64158, 64186, 64214, 64242, 64270, 64309, 64338, 64366, 64395, 64422, 64450, 64478, 64506, 64534, 64562, 64590, 64618, 64647, 64674, 64703, 64731, 64759, 64787, 64815, 64853, 64882, 64910, 64938, 64966, 64994, 65022, 65050, 65078, 65107, 65134, 65162, 65191, 65218, 65247, 65275, 65303, 65331, 65369, 65398, 65426, 65454, 65482, 65510, 65538, 65566, 65594, 65623, 65650, 65679, 65707, 65735, 65763, 65802, 65830, 65858, 65898, 65925, 65953, 65994, 66022, 66050, 66078, 66106, 66134, 66162, 66190, 66219, 66246, 66274, 66302, 66330, 66358, 66386, 66414, 66442, 66471, 66498, 66526, 66555, 66583, 66611, 66639, 66666, 66695, 66722, 66751]
                xData = [-165.9, -163.5, -161.4, -157.5, -153.6, -150.45, -146.4, -143.4, -140.25, -136.35, -134.1, -130.95, -126.75, -124.35, -121.95, -118.8, -116.4, -113.55, -111.45, -110.55, -108.6, -105.9, -104.4, -102.9, -101.7, -100.8, -99.15, -98.1, -97.65, -97.5, -96.6, -94.95, -94.2, -94.35, -93.9, -93.3, -92.85, -92.85, -94.35, -94.95, -95.1, -95.7, -96.3, -97.2, -98.25, -98.4, -99, -100.2, -102.45, -105.75, -108.9, -111.15, -115.2, -119.25, -123.3, -126.9, -130.05, -133.95, -138.45, -142.8, -147, -151.35, -152.25, -153.45, -153.15, -154.35, -154.2, -153.3, -153.6, -156.3, -158.1, -162.75, -164.85, -168, -172.5, -175.35, -179.1, -183.9, -186.45, -190.35, -193.5, -196.2, -197.85, -199.65, -201, -202.95, -204.6, -205.05, -205.95, -209.1, -211.05, -212.25, -214.35, -215.55, -219.3, -222.6, -225.75, -229.5, -234.6, -238.5, -243.45, -247.5, -251.4, -254.7, -257.25, -260.4, -262.05, -265.05, -266.4, -267.45, -270.15, -272.55, -274.2, -276, -277.35, -280.05, -283.95, -285.15, -286.5, -287.55, -288.9, -289.8, -290.1, -289.65, -289.95, -290.7, -291.3, -292.95, -294.45, -295.65, -297, -297.9, -298.95, -300, -300.6, -301.35, -301.65, -302.4, -303.6, -305.4, -305.55, -306.9, -307.65, -308.85, -308.1, -307.95, -307.35, -308.7, -307.65, -306.45, -306.45, -307.35, -307.2, -306.9, -305.4, -305.85, -306.75, -306.9, -306.15, -305.85, -306.15, -306.45, -306, -304.95, -303.75, -303.9, -303.9, -302.85, -301.95, -300.75, -300.15, -300, -298.5, -297.9, -296.85, -296.1, -295.5, -294.75, -292.65, -291.75, -290.25, -287.85, -285.45, -282.9, -280.2, -277.95, -274.8, -270, -267.75, -265.2, -262.2, -259.65, -257.1, -254.25, -251.25, -247.95, -244.8, -242.85, -238.8, -235.5, -233.1, -231.9, -229.8, -227.55, -224.7, -223.5, -221.1, -218.55, -215.7, -213.6, -211.8, -208.95, -205.65, -202.5, -199.5, -197.4, -192.9, -189.45, -186.6, -183.45, -180, -176.55, -173.4, -170.7, -167.25, -164.25, -162, -160.2, -157.35, -153.45, -150.6, -148.2, -145.2, -141.45, -139.2, -136.2, -134.4, -132.45, -129.9, -129.3, -128.7, -126.15, -126.3, -124.65, -123.9, -123.45, -121.8, -120.75, -121.65, -121.05, -121.5, -121.8, -120.9, -119.7, -119.85, -117.15, -115.8, -114.9, -112.95, -111.45, -110.25, -107.25, -106.05, -104.25, -102, -100.5, -99.3, -98.25, -96.9, -95.7, -94.8, -93.9, -93.15, -92.7, -90.9, -90.45, -90.9, -90, -89.7, -89.1, -88.95, -88.65, -88.2, -86.4, -85.95, -85.35, -84.6, -83.85, -84.45, -84.3, -83.7, -84, -83.7, -82.8, -81.9, -80.85, -80.4, -81.45, -80.85, -79.8, -81.9, -83.55, -83.85, -85.35, -85.35, -85.95, -87, -86.25, -85.65, -85.95, -85.8, -86.85, -87.6, -89.25, -90.15, -90.6, -92.1, -93.45, -95.25, -96, -97.95, -99.6, -101.25, -104.1, -106.65, -108.9, -112.05, -113.85, -117.45, -120.9, -123.6, -126.15, -127.95, -131.7, -135.45, -137.55, -141.6, -145.35, -149.1, -153.9, -157.05, -160.2, -164.25, -167.55, -169.65, -171.3, -173.25, -175.2, -177.6, -180.3, -182.85, -185.7, -189.6, -193.65, -198.6, -204.3, -209.7, -214.8, -221.7, -229.2, -236.4, -243.45, -249, -255.3, -261.75, -266.1, -269.85, -272.7, -275.25, -277.05, -278.55, -280.05, -282, -282.9, -284.85, -285.6, -286.5, -287.1, -288.15, -290.1, -291.9, -292.2, -293.1, -295.8, -297.75, -299.25, -301.5, -301.5, -303.3, -304.8, -305.1, -305.4, -305.7, -304.05, -306.75, -307.2, -307.95, -307.95, -308.55, -308.1, -308.7, -306.45, -306.3, -306, -306.3, -307.2, -307.65, -307.8, -308.55, -308.4, -309.3, -309, -308.4, -307.95, -307.5, -306.45, -306.45, -304.8, -304.35, -303.15, -302.1, -300.9, -300.75, -298.95, -297, -294.15, -292.35, -290.7, -289.05, -287.7, -285.45, -283.5, -283.2, -280.8, -279, -277.5, -274.8, -273.45, -272.7, -271.35, -270.9, -269.1, -266.7, -264.45, -262.05, -257.55, -253.05, -249.3, -246.45, -244.05, -241.2, -238.35, -237.15, -235.8, -233.1, -229.5, -226.8, -224.85, -221.85, -218.4, -214.8, -211.2, -209.25, -204.6, -199.05, -195.15, -189.9, -184.2, -179.25, -173.1, -168.6, -164.25, -159.15, -155.7, -152.1, -148.35, -145.5, -142.05, -139.05, -136.5, -133.5, -132.3, -131.1, -128.55, -126.45, -124.95, -124.35, -124.05, -122.4, -120.75, -120.45, -120.3, -120.45, -119.55, -118.95, -117.6, -116.7, -115.5, -114, -112.05, -110.1, -108.45, -106.95, -105.45, -103.2, -102.3, -99.6, -97.8, -96, -94.5, -92.55, -90.45, -88.5, -88.95, -87.45, -86.55, -86.4, -86.55, -87.45, -87.15]
                yData = [2.55, 2.55, 3.15, 4.5, 4.05, 4.35, 5.25, 4.8, 3.75, 4.2, 3.9, 4.2, 3.9, 3.6, 3.3, 3.9, 2.4, 2.25, 1.5, 0.6, -1.05, -2.85, -4.35, -5.4, -8.1, -8.7, -10.05, -10.65, -11.4, -13.05, -14.55, -15, -18.3, -20.7, -22.95, -25.2, -27.15, -30.45, -34.05, -36.6, -38.25, -40.8, -43.2, -47.1, -48.6, -51, -53.55, -57.15, -61.05, -64.8, -67.35, -71.85, -76.65, -81.45, -85.95, -90.45, -94.35, -98.25, -102.3, -105.9, -109.35, -112.05, -113.55, -115.2, -115.8, -115.2, -113.25, -112.2, -111.45, -111.75, -113.1, -114.45, -116.4, -120.3, -123.45, -126.75, -129.6, -131.55, -133.95, -135.6, -136.5, -136.95, -138.45, -139.5, -140.25, -141.9, -143.25, -144.15, -144.6, -144.6, -144.75, -145.65, -145.8, -146.1, -146.4, -147.9, -148.35, -148.8, -148.95, -148.8, -149.85, -150.75, -151.5, -152.4, -152.55, -152.4, -153.3, -153, -153.15, -151.65, -150.9, -152.25, -152.4, -152.25, -152.1, -151.35, -151.35, -151.2, -149.1, -148.2, -146.55, -145.8, -145.35, -146.1, -144.9, -144.9, -144.75, -144.45, -144, -142.95, -140.85, -141, -139.65, -139.35, -139.05, -137.7, -136.65, -136.35, -134.25, -133.95, -131.55, -130.05, -129.9, -129.3, -127.5, -127.05, -125.55, -125.25, -124.5, -123, -121.95, -120.75, -119.4, -119.25, -118.35, -117.75, -115.65, -113.85, -111.6, -110.7, -107.85, -106.35, -103.2, -101.55, -99.3, -98.25, -95.7, -94.35, -91.65, -91.05, -89.4, -89.1, -87.9, -87.45, -85.65, -85.2, -83.55, -82.2, -80.25, -78.15, -76.2, -74.4, -70.65, -68.1, -66, -61.95, -59.4, -54.9, -51.75, -50.1, -47.1, -44.1, -42.6, -41.1, -39.45, -37.2, -34.8, -32.7, -30.9, -28.65, -25.2, -23.4, -22.2, -21.45, -20.1, -18, -17.1, -17.4, -16.5, -15.3, -13.8, -13.05, -12.6, -10.8, -8.85, -7.65, -5.1, -2.85, -0.6, 0.75, 2.25, 3.6, 4.95, 5.7, 6.15, 7.05, 7.65, 8.25, 9, 9.6, 10.05, 11.4, 11.7, 13.05, 13.35, 15, 15.6, 15.75, 14.7, 14.55, 14.85, 14.7, 13.65, 13.65, 13.8, 14.4, 14.55, 14.1, 14.7, 14.55, 14.55, 15, 15, 14.85, 14.85, 14.1, 13.8, 13.5, 12.3, 11.55, 12.3, 12, 12.45, 12, 11.1, 10.8, 11.1, 10.65, 10.35, 10.65, 10.65, 10.95, 9.9, 9.3, 7.65, 6.6, 5.25, 4.2, 2.4, 1.95, 1.35, 0.9, 0.45, -1.5, -2.85, -3.45, -3.15, -4.95, -6.75, -8.7, -10.2, -11.4, -12.75, -15.6, -16.05, -16.65, -16.8, -17.55, -18.45, -19.65, -19.8, -21.45, -22.65, -25.05, -26.85, -28.35, -30.45, -32.4, -33.15, -35.55, -37.05, -38.25, -39.6, -39.6, -41.7, -45.3, -46.8, -48.9, -51.15, -54.15, -58.35, -61.05, -63.3, -66.3, -69, -71.7, -75, -77.85, -81, -84.6, -87.6, -91.5, -95.55, -97.8, -100.95, -103.8, -106.5, -108.75, -110.85, -112.65, -115.35, -117, -119.4, -121.65, -124.2, -126, -127.5, -129.45, -131.55, -132.9, -133.8, -135, -136.95, -139.2, -141.3, -142.8, -144.45, -147.15, -149.85, -151.05, -152.25, -152.85, -153.3, -154.8, -154.05, -154.05, -154.8, -154.65, -153.75, -153.15, -152.4, -152.55, -151.95, -150.6, -150.3, -149.55, -149.7, -149.4, -148.2, -147.3, -147.3, -146.55, -146.25, -144.75, -143.1, -142.05, -140.55, -138, -135.9, -135, -133.95, -133.05, -131.1, -130.8, -130.35, -129.9, -129.3, -128.4, -127.35, -127.35, -125.7, -124.65, -123.75, -122.55, -122.7, -121.65, -120.75, -119.7, -119.1, -118.35, -116.85, -114.15, -112.95, -111.3, -110.4, -108.75, -106.8, -105, -103.65, -101.4, -98.4, -94.8, -91.05, -87.9, -84.6, -81.75, -78.3, -76.35, -74.7, -73.8, -71.7, -69.75, -68.55, -67.35, -65.85, -63.6, -62.7, -61.95, -61.2, -59.1, -57.9, -56.4, -54.45, -51.3, -48.6, -45.3, -42.6, -39.3, -35.85, -34.35, -32.55, -31.65, -30.75, -28.65, -28.2, -27.3, -25.65, -23.55, -20.85, -19.8, -18.45, -15.75, -14.4, -12, -9, -6.15, -2.4, -0.15, 1.65, 3.75, 5.55, 5.55, 6.45, 7.05, 8.4, 9.75, 10.95, 11.7, 13.8, 14.25, 14.55, 13.95, 14.25, 14.1, 13.95, 13.5, 14.4, 14.1, 14.4, 13.35, 13.35, 13.05, 12.15, 10.95, 11.25, 11.7, 11.7, 12.15, 11.85, 12.75, 12.9, 13.05, 12.9, 12.45, 11.4, 10.8, 10.5, 8.85, 7.5, 6.15, 4.95, 4.35, 3.15, 1.05, 0.6, -1.05, -2.7, -3, -5.25, -6.3]
                zData = [379.2, 377.7, 374.7, 369.15, 366.15, 362.85, 359.25, 356.85, 353.4, 349.8, 346.8, 343.65, 338.85, 334.8, 330.6, 325.5, 321, 316.8, 312.9, 310.2, 307.95, 303, 300.9, 297.45, 294.9, 291.6, 289.65, 286.65, 285.15, 282.3, 280.8, 277.95, 274.8, 269.85, 265.5, 262.5, 258.6, 254.4, 252.45, 249.6, 247.95, 247.05, 245.85, 244.65, 243.9, 242.55, 242.1, 241.35, 240.45, 239.55, 237.45, 236.55, 234.6, 234, 233.7, 233.85, 233.7, 235.65, 236.1, 238.35, 238.35, 240.15, 240.75, 240.6, 240.75, 241.65, 240.9, 241.2, 241.2, 241.8, 244.2, 245.55, 247.65, 249.9, 253.65, 255.15, 258.45, 260.25, 262.65, 263.7, 265.8, 266.25, 267.9, 267.9, 269.1, 270.45, 272.55, 273, 274.95, 275.7, 277.8, 279.3, 281.4, 283.65, 287.4, 290.1, 294.9, 298.2, 302.25, 305.85, 309, 312.6, 316.2, 319.35, 322.65, 326.4, 328.8, 332.25, 333.9, 337.05, 338.25, 339.75, 340.35, 342, 343.5, 345.75, 346.95, 350.1, 352.8, 356.1, 358.05, 359.55, 362.25, 364.65, 366.3, 368.7, 369.45, 371.85, 373.8, 374.85, 375.75, 377.85, 378.45, 380.25, 382.8, 385.35, 387.45, 388.65, 390.3, 391.2, 392.4, 392.55, 393.75, 394.65, 397.5, 397.95, 400.95, 402.15, 403.8, 403.5, 405.3, 406.35, 408.3, 409.2, 410.85, 412.05, 415.35, 416.85, 418.05, 419.1, 421.5, 423.75, 426.3, 427.05, 428.1, 429.75, 431.4, 431.1, 432.3, 431.4, 433.05, 433.8, 435.15, 435, 435.3, 435.45, 437.25, 437.1, 437.1, 436.35, 436.8, 437.7, 437.7, 437.7, 437.85, 438.3, 438, 438.15, 437.4, 436.65, 436.05, 435.3, 433.65, 433.65, 432.6, 432.15, 431.25, 430.35, 429.6, 429.6, 427.05, 425.85, 424.35, 422.7, 421.5, 420.6, 418.65, 418.8, 417.45, 416.55, 416.55, 414.45, 412.8, 411.3, 408.3, 406.35, 403.8, 400.95, 398.55, 396, 393.6, 390.6, 388.8, 386.7, 383.7, 381.6, 378.3, 376.2, 374.85, 371.4, 368.55, 365.7, 363, 360.6, 357.75, 355.5, 353.25, 351.75, 350.55, 349.05, 348.75, 347.25, 345.9, 345.45, 345, 345.15, 344.1, 342.15, 342.45, 342.15, 342, 340.35, 338.55, 337.8, 337.2, 335.25, 333.75, 331.65, 329.25, 327, 324.75, 321.75, 319.35, 317.55, 314.85, 312.6, 310.5, 308.4, 306.6, 305.25, 302.1, 301.5, 300.9, 299.55, 298.35, 297.3, 296.1, 294.9, 291.45, 288.9, 286.5, 285.15, 283.8, 280.65, 278.7, 277.5, 276.6, 274.35, 272.7, 269.85, 268.8, 267.45, 265.65, 263.55, 262.8, 260.7, 260.85, 259.65, 258.15, 257.4, 256.35, 254.25, 252.45, 249.75, 249.3, 248.4, 247.5, 246.45, 245.85, 244.8, 244.2, 241.95, 240.6, 239.7, 238.2, 237.9, 237, 235.5, 234.45, 233.25, 232.2, 232.05, 231.45, 231.6, 232.2, 231.45, 231.75, 231.45, 229.95, 230.4, 228.9, 228.9, 231, 232.35, 232.2, 234.75, 235.65, 238.5, 238.95, 240.3, 241.35, 244.05, 246.45, 248.4, 249.6, 252.3, 252.6, 253.65, 254.85, 254.7, 256.35, 259.05, 261.6, 266.1, 269.85, 274.95, 279.9, 284.7, 289.5, 295.35, 299.7, 305.4, 310.2, 316.05, 322.5, 327.75, 333, 337.65, 340.8, 344.7, 347.4, 349.05, 351.45, 352.95, 354.3, 356.7, 358.2, 360.75, 361.95, 363.9, 365.55, 368.25, 370.05, 372.75, 375.45, 379.05, 382.5, 385.05, 387.3, 390.45, 393.15, 394.35, 396.45, 396.9, 398.4, 400.35, 401.85, 401.85, 403.65, 403.65, 404.1, 404.7, 404.55, 405.45, 406.5, 406.95, 407.55, 409.35, 410.55, 412.8, 413.25, 414.15, 416.1, 416.7, 417.45, 418.95, 420, 421.2, 423.3, 424.35, 427.8, 430.35, 431.7, 434.25, 436.8, 438.45, 439.2, 440.1, 440.4, 441.45, 440.25, 439.95, 439.5, 440.1, 440.25, 440.7, 440.85, 442.35, 442.35, 442.35, 442.5, 441.6, 441.3, 439.8, 438.6, 437.4, 436.2, 434.1, 433.05, 430.5, 429.9, 428.25, 427.05, 426.9, 425.7, 423.9, 423.45, 421.05, 418.5, 416.1, 413.7, 410.85, 408.45, 405.6, 403.05, 400.65, 398.25, 394.65, 391.8, 388.5, 384.6, 381.6, 378.6, 375.45, 371.1, 368.4, 366, 363.6, 361.35, 359.25, 356.7, 354.6, 353.1, 350.85, 350.1, 348.75, 347.4, 346.5, 347.1, 345.45, 344.55, 343.05, 340.65, 339, 337.2, 335.7, 334.05, 332.25, 328.2, 325.95, 323.55, 321.6, 318, 314.55, 311.85, 310.2, 307.95, 305.55, 302.25, 300.3, 298.5, 295.8, 294.3, 291.9, 288.6, 286.05, 283.2, 281.1, 277.35, 274.35]
                uData = [-829.05, -829.80, -828.00, -837.45, -836.40, -856.65, -862.20, -860.85, -860.10, -868.95, -862.50, -861.15]
                vData = [-367.35, -360.60, -181.65, -435.30, -251.85, -283.35, -288.00, -284.40, -317.25, -341.40, -237.30, -209.55]
                break
        }
        
        // transpose the three arrays into array of triples
        for (let n = 0; n < scanTimes.length; n++){
            let xyz = []
            xyz.push(xData[n])
            xyz.push(yData[n])
            xyz.push(zData[n])
            scanData.push(xyz)
        }
    }

}
