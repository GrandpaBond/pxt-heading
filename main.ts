/**
 * An extension providing a compass-heading for a buggy located anywhere on the globe 
 * (except at the magnetic poles!), with any mounting orientation for its microbit.
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
    //const Window = 100 // minimum spacing (in samples) of Ellipse major-axis candidates

    // CLASSES

    // Characteristics of the Ellipse formed when projecting the Spin-Circle onto a View plane
    class Ellipse {
        plane: string; // name (for debug)
        uDim: number; // horizontal axis of this View
        vDim: number; // vertical axis of this View
        uOff: number; // horizontal offset needed to re-centre the Ellipse along the U-axis
        vOff: number; // vertical offset needed to re-centre the Ellipse along the V-axis
        // properties used in the Spin-Circle scan:
        rSq: number;  // latest sample's smoothed radius-squared
        angle: number; // latest sample's polar angle on the Ellipse 
        angleDegrees: number; // testing
        hiRsq: number; // max radius-squared 
        loRsq: number; // min radius-squared
        // history:
        angleWas: number; // previous sample's angle
        rSqWas: number; // previous sample's rSq
        rSqWasWas: number // last-but-one sample's rSq
        // projection params from latest detected major axis:
        theta: number; // angle (in radians) anticlockwise from U-axis to major axis 
        thetaDegrees: number; // testing
        cosTheta: number; // save for efficiency
        sinTheta: number; // ditto
        scale: number; // stretch in V needed to make Ellipse circular again
        // others:
        fromBelow: number; // rotation seen by this view of the clockwise scan
        firstAngle: number; // angle of first scan sample 
        turned: number; // scan angle accumulator (in radians)
        period: number; // this View's assessment of average rotation time


        constructor(plane: string, dim0: number, dim1: number, off0: number, off1: number) {
            this.plane = plane // (as a DEBUG aid)
            this.uDim = dim0
            this.vDim = dim1
            this.uOff = off0
            this.vOff = off1
            //this.peaks = []       
        }

        // convert the raw values to (kind of) re-centred polar coordinates (rSq,angle)
        // (no need to extract sqrt of rSq every time)
        polarise(uRaw: number, vRaw: number) {
            let u = uRaw - this.uOff
            let v = vRaw - this.vOff
            this.rSq = (u * u) + (v * v)
            this.angle = Math.atan2(v,u)
            this.angleDegrees = radians2degrees(this.angle)
        }
  
        // initialise the smoothing and turning history
        firstSample(uRaw: number, vRaw:number) {
            this.polarise(uRaw, vRaw) // derive (rSq,angle) for first sample
            this.angleWas = this.angle
            this.firstAngle = this.angle
            this.turned = 0
            this.rSqWasWas = this.rSq
            this.rSqWas = this.rSq
            this.hiRsq = this.rSq
            this.loRsq = this.rSq
        } 

        // process another scan sample
        nextSample(index: number, uRaw: number, vRaw: number) {
            this.polarise(uRaw, vRaw) // derive (rSq, angle) for this sample
            // While tracking the radius, we use inertial smoothing to reduce multiple
            // spurious peak-detections due to minor fluctuations in readings
            let ratio = this.rSqWas / this.rSqWasWas
            let guess = this.rSqWas * ratio // projecting same grow/shrink factor
            this.rSq = (guess * Inertia) + (this.rSq * (1 - Inertia)) // deflect by new value

            // If slope changes from rising to falling we've potentially just passed a major-axis.
            if (   (ratio >= 1) // was rising or level...
                && (this.rSq < this.rSqWas)) { // ... but now definitely falling
                // New peak! if radius is longest so far, adopt it as definitive for now
                if (this.rSqWas > this.hiRsq) {
                    this.hiRsq = this.rSqWas
                    this.theta = this.angleWas
                    this.thetaDegrees = radians2degrees(this.theta) // (while testing)
 
                }
            }
            // also keep track of the smallest minor-axis radius
            this.loRsq = Math.min(this.loRsq, this.rSq) // shortest so far...

            // keep track of complete revolutions, by detecting roll-round
            if ((this.angleWas < 0) && (this.angle > 0)) this.turned += 2 * Math.PI
            this.angleWas = this.angle   

            if (logging && testing) {
                datalogger.log(
                    datalogger.createCV("view", this.plane),
                    datalogger.createCV("index", index),
                    datalogger.createCV("loRsq", this.loRsq),
                    datalogger.createCV("rSq", this.rSq),
                    datalogger.createCV("hiRsq", this.hiRsq),
                    datalogger.createCV("thetaDegrees", this.thetaDegrees),
                    datalogger.createCV("turned", this.turned))
            }

            // update smoothing history
            this.rSqWasWas = this.rSqWas
            this.rSqWas = this.rSq
        }

        // perform final analysis after scanning
        finish() {        
            // Preset the rotation factors for aligning major axis with the U-axis
            this.cosTheta = Math.cos(this.theta)
            this.sinTheta = Math.sin(this.theta)
            // Use the ratio of the detected extreme radii to derive the V-axis scaling factor
            // needed to stretch the Ellipse back into a circle            
            this.scale = Math.sqrt(this.hiRsq / this.loRsq)
            // Work out the total scanned angle and derive the rotation period
            let revs = (this.angle + this.turned - this.firstAngle) / (2 * Math.PI)
            this.period = scanTime / revs
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

            return angle
        }
    }

    // GLOBALS
   
    let scanTimes: number[]  = [] // sequence of time-stamps for scanned readings 
    let scanData: number[][] = [] // scanned sequence of [X,Y,Z] magnetometer readings
    let scanTime: number = 0 // duration of scan in ms
    let views: Ellipse[] = [] // the three possible elliptical views of the Spin-Circle
    let bestView = -1
    let uDim = -1 // the best "horizontal" axis (pointing East) for transformed readings (called U)
    let vDim = -1 // the best "vertical" axis (pointing North) for transformed readings (called V)
    let toNorth = 0 // the angular bias to be added (so that North = 0)
    let toNorthDegrees = 0 // (testing)
    let strength = 0 // the average magnetic field-strength observed by the magnetometer
    let fromBelow = false // set "true" if orientation means readings project backwards

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
        }
        if (testing) {
           simulateScan("Y-up") // Y-Axis vertical; spun around Y-Axis
            // simulateScan("tetra")
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
        
        if (logging && !testing) {

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

        // we'll typically need at least a couple of second's worth of scanned readings...
        let nSamples = scanTimes.length
        scanTime = scanTimes[nSamples - 1] - scanTimes[0]
        if (scanTime < 2000) {
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
            xlo = Math.min(xlo, scanData[i][Dim.X])
            yhi = Math.max(yhi, scanData[i][Dim.Y])
            ylo = Math.min(ylo, scanData[i][Dim.Y])
            zhi = Math.max(zhi, scanData[i][Dim.Z])
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


        // process scan from 2nd sample onwards...
        for (let j = 1; j < nSamples; j++) {
            xRaw = scanData[j][Dim.X]
            yRaw = scanData[j][Dim.Y]
            zRaw = scanData[j][Dim.Z]
            views[View.XY].nextSample(j, xRaw, yRaw) // projection in XY plane
            views[View.YZ].nextSample(j, yRaw, zRaw) // projection in YZ plane
            views[View.ZX].nextSample(j, zRaw, xRaw) // projection in ZX plane

            // Square of overall 3-D field strength is sum of squares in each View
            // Here, we will end up accumulating each one twice over!
            strength += views[View.XY].rSq
            strength += views[View.YZ].rSq
            strength += views[View.ZX].rSq
        }

         // check average overall field-strength (undoing the double-counting)
        strength = Math.sqrt((strength / 2) / nSamples)
        if (strength < MarginalField) {
            return -2  // "FIELD STRENGTH TOO WEAK"
        }

        // for each view, set up projection metrics, using the most recently detected major axis
        // Also derive rotation period (should match between views!)
        views[View.XY].finish()
        views[View.YZ].finish()
        views[View.ZX].finish()
   
       // check that at least one View saw a complete rotation (= 2*Pi radians)...
        if (   (views[View.XY].turned < (2 * Math.PI))
            && (views[View.YZ].turned < (2 * Math.PI))
            && (views[View.ZX].turned < (2 * Math.PI)) ) {
            return -3 // "NOT ENOUGH SCAN ROTATION"
        }

        // Choose the "roundest" Ellipse, the one with lowest (eccentricity = scale).
        bestView = View.XY
        if (views[View.YZ].scale < views[bestView].scale) bestView = View.YZ
        if (views[View.ZX].scale < views[bestView].scale) bestView = View.ZX

        basic.clearScreen()
        basic.pause(100)
        basic.showString(views[bestView].plane)
        basic.pause(300)
        
        // Depending on mounting orientation, the bestView Ellipse might possibly be seeing the 
        // Spin-Circle from "underneath", effectively resulting in an anti-clockwise scan.
        // Check polarity of turns:
        fromBelow = views[bestView].turned > 0

        // extract some bestView fields into globals for future brevity and efficiency...
        uDim = views[bestView].uDim
        vDim = views[bestView].vDim
        
        // we've now finished with the scanning data, so release its memory
        scanTimes = []
        scanData = []

        // SUCCESS!
        return 0
    }



    /**
     * Read the magnetometer and register the buggy's current direction as "North",
     * (i.e. the direction that will in future return zero as its heading).
     * The actual direction the buggy is pointing when this function is called could be
     * Magnetic North; True North (compensating for local declination); or any convenient
     * direction from which to measure subsequent heading angles.
     */
    //% block="set North" 
    //% inlineInputMode=inline 
    //% weight=80 
    export function setNorth() {
        // We have successfully set up the projection parameters. Now we need to relate them to "North".
        // Take the average of seven new readings to get a stable fix on the current heading
        let uRaw = 0
        let vRaw = 0
        if (testing) { //arbitrarily choose first test reading for North
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

        // get the projection angle of the [uRaw,vRaw] vector on the Spin-Circle, 
        // and remember this as the (global) fixed bias to North
        toNorth = views[bestView].project(uRaw, vRaw)
        toNorthDegrees = radians2degrees(toNorth)


        if (logging) {
            datalogger.setColumnTitles("uDim", "vDim", "uOff", "vOff",
              "thetaDegrees", "scale", "period", "toNorthDegrees", "strength")

            datalogger.log(
                datalogger.createCV("uDim", uDim),
                datalogger.createCV("vDim", vDim),
                datalogger.createCV("uOff", views[bestView].uOff),
                datalogger.createCV("vOff", views[bestView].vOff),
                datalogger.createCV("theta", views[bestView].thetaDegrees),
                datalogger.createCV("scale", views[bestView].scale),
                datalogger.createCV("period", views[bestView].period),
                datalogger.createCV("toNorthDegrees", toNorthDegrees),
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
            if (test > 7) test = 0 // roll round 8 points of the compass
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

        // project the reading from Ellipse-view to Spin-Circle, as radians anticlockwise from U-axis
        let onCircle = views[bestView].project(uRaw, vRaw)
        let onCircleDegrees = radians2degrees(onCircle)
        // subtract toNorth to get radians anticlockwise from North,
        // which is just what we wanted if viewed "fromBelow"
        let angle = onCircle - toNorth
        let angleDegrees = radians2degrees(angle)
        // otherwise, negate heading angle to measure clockwise from North.
        if (!fromBelow) angle = -angle
        if (logging) {
            datalogger.log(datalogger.createCV("uRaw", uRaw),
                datalogger.createCV("vRaw", vRaw),
                datalogger.createCV("onCircleDegrees", onCircleDegrees),
                datalogger.createCV("toNorthDegrees", toNorthDegrees),
                datalogger.createCV("angleDegrees", angleDegrees))
        }
        return radians2degrees(angle)
    }

    /**
     * Return average RPM of the most recent scan 
     */
    //% block="scanRPM" 
    //% inlineInputMode=inline 
    //% weight=50 
    export function scanRPM(): number {
        if (views[bestView].period <= 0) {
            return -4 // ERROR: SUCCESSFUL SCAN IS NEEDED FIRST
        } else {
            return 60000 / views[bestView].period
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
                scanTimes = [12068, 12096, 12124, 12152, 12180, 12208, 12236, 12264, 12292, 12322, 12348, 12376, 12404, 12432, 12460, 12488, 12516, 12544, 12572, 12600, 12628, 12656, 12684, 12712, 12740, 12768, 12796, 12824, 12852, 12883, 12908, 12936, 12964, 12992, 13020, 13048, 13076, 13104, 13132, 13160, 13188, 13216, 13244, 13272, 13300, 13328, 13359, 13384, 13412, 13440, 13468, 13496, 13524, 13552, 13580, 13608, 13636, 13664, 13692, 13720, 13748, 13776, 13807, 13832, 13860, 13888, 13916, 13944, 13972, 14000, 14028, 14056, 14084, 14112, 14140, 14171, 14196, 14224, 14252, 14280, 14308, 14336, 14364, 14392, 14420, 14448, 14476, 14504, 14535, 14560, 14588, 14616, 14644, 14672, 14701, 14728, 14756, 14784, 14812, 14844, 14872, 14900, 14928, 14956, 14984, 15012, 15041, 15072, 15100, 15128, 15156, 15184, 15212, 15240, 15268, 15297, 15328, 15356, 15384, 15412, 15440, 15468, 15496, 15524, 15556, 15584, 15612, 15640, 15668, 15696, 15724, 15752, 15784, 15812, 15840, 15868, 15896, 15924, 15952, 15984, 16012, 16040, 16068, 16097, 16124, 16156, 16184, 16213, 16240, 16268, 16300, 16329, 16356, 16384, 16412, 16444, 16472, 16500, 16528, 16556, 16588, 16616, 16644, 16672, 16704, 16732, 16760, 16788, 16820, 16848, 16876, 16908]
                xData = [-160.80, -184.05, -191.85, -198.90, -207.30, -213.60, -221.70, -227.85, -232.35, -237.75, -243.00, -248.25, -252.45, -255.90, -260.85, -265.05, -268.65, -272.70, -275.10, -279.00, -281.70, -283.20, -285.75, -287.55, -288.45, -289.95, -290.70, -291.15, -291.75, -291.15, -291.00, -289.95, -289.35, -288.00, -287.55, -286.05, -285.00, -283.05, -282.00, -279.90, -277.65, -274.50, -270.30, -266.40, -261.75, -256.05, -251.25, -247.05, -241.95, -237.30, -232.35, -227.40, -223.95, -218.55, -212.85, -207.45, -203.70, -197.85, -192.30, -186.90, -181.65, -175.65, -170.25, -163.05, -156.00, -150.60, -144.30, -138.15, -132.30, -126.45, -121.50, -116.55, -111.15, -105.15, -100.20, -95.85, -90.75, -84.60, -80.10, -75.30, -72.00, -68.10, -63.60, -59.55, -57.15, -54.90, -52.95, -50.55, -48.30, -47.40, -47.25, -47.85, -48.15, -48.90, -50.10, -53.10, -55.80, -58.95, -61.20, -64.35, -67.95, -72.15, -76.05, -80.40, -83.40, -87.75, -91.65, -94.20, -97.65, -102.00, -104.85, -109.80, -112.65, -117.90, -124.95, -131.10, -135.90, -142.35, -147.30, -153.90, -159.75, -164.70, -169.95, -174.60, -180.00, -186.75, -192.90, -197.85, -202.80, -206.40, -211.35, -214.80, -219.30, -223.20, -226.35, -229.80, -234.45, -238.20, -242.10, -244.80, -248.40, -253.20, -256.80, -260.55, -264.30, -267.30, -271.80, -274.65, -277.20, -280.50, -282.45, -285.45, -288.60, -289.50, -290.55, -291.45, -292.80, -292.35, -292.05, -292.05, -291.45, -291.45, -290.55, -287.70, -287.70, -286.80, -284.85, -281.70, -279.30, -275.85, -273.75, -269.55]
                yData = [-16.50, -16.65, -15.75, -14.55, -14.70, -14.25, -14.55, -14.25, -14.40, -14.70, -14.40, -14.25, -14.10, -13.95, -14.25, -14.10, -13.95, -14.40, -14.25, -15.15, -15.00, -15.00, -15.45, -15.90, -15.75, -16.20, -15.45, -15.45, -16.05, -15.15, -14.40, -14.55, -14.55, -14.70, -14.70, -13.80, -14.25, -15.00, -15.00, -14.70, -14.70, -14.55, -14.40, -14.70, -14.70, -14.40, -15.00, -14.85, -15.45, -15.45, -15.60, -15.15, -14.55, -13.95, -13.95, -13.50, -14.55, -15.15, -15.90, -16.95, -17.25, -17.70, -18.30, -17.70, -18.00, -17.10, -16.50, -17.10, -16.95, -16.05, -15.75, -14.40, -14.40, -14.85, -13.95, -14.10, -14.25, -15.00, -15.45, -15.15, -15.00, -15.00, -13.65, -13.95, -13.20, -12.75, -13.35, -13.35, -13.65, -14.40, -14.25, -15.30, -16.05, -16.35, -16.35, -16.65, -16.80, -16.80, -16.05, -16.05, -15.75, -15.60, -15.15, -15.30, -15.45, -15.75, -15.45, -14.85, -15.45, -15.60, -15.30, -14.85, -15.00, -16.05, -17.40, -16.50, -16.65, -17.40, -18.45, -18.45, -17.55, -16.65, -17.70, -17.25, -16.80, -16.20, -15.60, -15.30, -14.70, -13.35, -13.80, -14.55, -15.45, -14.55, -14.40, -14.55, -14.70 - 27.15, -26.85, -27.00, -27.75, -27.90, -27.00, -27.15, -26.55, -25.50, -24.60, -24.15, -25.05, -25.35, -24.90, -25.35, -26.25, -26.55, -26.70, -25.50, -25.95, -26.25, -25.50, -24.45, -25.05, -25.35, -26.25, -25.65, -25.95, -26.55, -27.60, -27.00, -27.30, -26.70, -26.85, -25.95, -25.80, -25.50, -27.15, -26.70, -26.55, -26.55, -27.00, -26.70, -26.70, -25.50, -25.65, -25.50, -25.50, -25.80, -26.70, -27.30, -27.30, -28.05, -28.20, -28.80, -29.55, -28.65, -28.50, -28.65, -28.80, -29.10, -28.95, -28.80, -30.60, -30.30, -31.05, -30.90, -31.05, -31.35, -31.35, -31.20, -31.65, -31.50, -31.65, -31.65, -31.65, -31.80, -30.60, -30.45, -30.30, -29.70, -30.45, -30.30, -29.55, -29.70, -30.15, -30.45, -30.30, -29.85, -30.00, -30.30, -30.90, -30.60, -30.45, -31.05, -31.35, -31.05, -30.90, -31.05, -30.75, -29.85, -29.85, -30.15, -30.75, -30.15, -28.65, -28.80, -29.10, -29.25, -29.25, -28.05, -28.80, -29.85, -29.40, -29.10, -28.95, -28.20, -28.35, -27.75, -27.75, -28.35, -28.50, -28.80, -28.65, -29.10, -28.95, -29.25, -28.50, -28.65, -27.60, -27.45, -27.00, -27.15, -26.10, -26.10, -25.95, -26.10, -26.55, -26.70, -27.75, -28.50, -28.35, -28.50, -28.65, -28.35, -28.65, -27.30, -26.70, -27.00, -26.70, -26.25, -26.10, -25.35, -25.65, -25.05, -25.65, -25.80, -25.35, -25.35, -26.25, -27.00, -27.30, -26.70, -26.85, -27.60, -27.75, -26.85, -26.70, -26.85, -26.85, -27.75, -28.20]
                zData = [-168.90, -191.25, -190.20, -189.00, -187.20, -184.20, -180.45, -178.65, -176.70, -173.55, -170.55, -166.50, -163.65, -160.95, -155.85, -149.40, -144.60, -138.90, -134.70, -127.50, -120.75, -115.95, -111.75, -105.90, -99.30, -93.30, -88.65, -84.60, -78.45, -73.20, -67.35, -63.30, -57.60, -51.60, -45.30, -40.80, -35.40, -30.60, -26.25, -22.05, -19.50, -15.45, -9.90, -4.65, -0.90, 4.65, 10.20, 15.90, 19.95, 22.95, 25.80, 28.95, 31.20, 32.55, 34.35, 37.65, 39.90, 41.85, 43.05, 44.25, 45.90, 46.80, 45.75, 45.60, 44.55, 44.85, 44.25, 42.15, 40.20, 38.85, 37.65, 35.10, 31.20, 28.80, 26.25, 23.55, 20.10, 15.15, 11.10, 7.95, 1.35, -4.35, -10.50, -16.35, -22.35, -29.10, -36.15, -42.90, -49.95, -58.35, -66.60, -74.55, -82.35, -90.90, -98.55, -105.90, -112.35, -118.20, -125.10, -130.65, -135.90, -139.20, -144.15, -148.65, -153.75, -158.10, -161.10, -164.10, -168.00, -170.55, -172.35, -174.30, -175.35, -178.50, -180.60, -183.75, -186.75, -188.70, -189.15, -191.10, -190.95, -191.10, -190.80, -190.05, -189.90, -190.05, -188.40, -187.95, -188.10, -187.80, -186.60, -184.65, -184.50, -183.15, -181.50, -178.05, -174.75, -172.35, -169.50, -165.45, -161.85, -158.70, -155.55, -152.10, -147.00, -143.70, -139.65, -134.55, -129.60, -124.65, -119.85, -114.75, -109.35, -104.40, -100.50, -93.90, -87.90, -81.00, -76.20, -69.75, -63.90, -58.05, -54.00, -49.65, -45.45, -40.05, -35.55, -30.75, -26.25, -21.30, -16.80, -11.40]
                uData = [-14.85, -104.1, -174.9, -188.55, -126.15, -39.6, 33.00, 42.15, -13.8]
                vData = [-57.30, -51.45, -113.55, -205.35, -277.2, -283.65, -226.50, -127.80, -61.2]
                break
            default: // tetrahedral: spin-axis 45 degrees to X, Y & Z
                scanTimes = [51901, 51929, 51957, 51985, 52013, 52041, 52069, 52097, 52125, 52155, 52182, 52209, 52237, 52265, 52293, 52321, 52349, 52377, 52405, 52433, 52461, 52489, 52517, 52545, 52573, 52601, 52629, 52657, 52685, 52716, 52741, 52769, 52797, 52825, 52853, 52881, 52909, 52937, 52965, 52993, 53021, 53049, 53077, 53105, 53134, 53161, 53192, 53217, 53245, 53273, 53301, 53329, 53357, 53385, 53413, 53441, 53469, 53497, 53525, 53553, 53581, 53609, 53640, 53669, 53697, 53725, 53754, 53781, 53809, 53837, 53865, 53893, 53921, 53950, 53978, 54009, 54037, 54065, 54093, 54121, 54149, 54177, 54205, 54233, 54262, 54289, 54318, 54346, 54377, 54405, 54433, 54461, 54489, 54517, 54545, 54573, 54601, 54630, 54658, 54686, 54716, 54741, 54769, 54797, 54825, 54853, 54882, 54913, 54941, 54969, 54997, 55025, 55054, 55081, 55110, 55138, 55166, 55196, 55221, 55249, 55277, 55306, 55334, 55362, 55390, 55421, 55449, 55477, 55505, 55533, 55562, 55590, 55621, 55649, 55677, 55705, 55734, 55762, 55790, 55821, 55849, 55877, 55906, 55934, 55962, 55993, 56021, 56049, 56077, 56106, 56137, 56165, 56193, 56221, 56250, 56278, 56309, 56337, 56366, 56394, 56425, 56453, 56482, 56510, 56541, 56569, 56598, 56626, 56657, 56686, 56714, 56742, 56773, 56802, 56830, 56861, 56890, 56921, 56949, 56977, 57005, 57038, 57065, 57094, 57122, 57150, 57178, 57206, 57234, 57262, 57290, 57318, 57346, 57374, 57402, 57430, 57458, 57486, 57514, 57542, 57570, 57598, 57626, 57654, 57682, 57710, 57744, 57770, 57798, 57826, 57854, 57882, 57910, 57938, 57966, 57994, 58022, 58050, 58078, 58106, 58134, 58162, 58190, 58218, 58246, 58274, 58308, 58333, 58362, 58390, 58418, 58446, 58474, 58502, 58530, 58558, 58586, 58614, 58642, 58670, 58698, 58726, 58754, 58782, 58810, 58845, 58874, 58902, 58930, 58958, 58986, 59014, 59042, 59070, 59098, 59126, 59154, 59182, 59210, 59238, 59266, 59294, 59329, 59358, 59386, 59414, 59442, 59470, 59498, 59526, 59554, 59582, 59610, 59638, 59666, 59694, 59722, 59757, 59786, 59814, 59842, 59870, 59898, 59926, 59954, 59982, 60010, 60038, 60066, 60094, 60122, 60150, 60185, 60214, 60249, 60278, 60306, 60334, 60362, 60390, 60418, 60446, 60474, 60502, 60530, 60558, 60593, 60621, 60650, 60678, 60706, 60734, 60762, 60790, 60818, 60853, 60881, 60910, 60938, 60966, 60994, 61022, 61050, 61078, 61113, 61142, 61170, 61198, 61226, 61254, 61282, 61310, 61345, 61374, 61402, 61430, 61458, 61486, 61514, 61549, 61578, 61606, 61634, 61662, 61690, 61725, 61754, 61782, 61810, 61838, 61873, 61902, 61930, 61958, 61986, 62021, 62050, 62078, 62106, 62134, 62169, 62198, 62226, 62254, 62289, 62318, 62346, 62375, 62409, 62438, 62466, 62502, 62530, 62558, 62594, 62622, 62650, 62685, 62713, 62741, 62769, 62797, 62825, 62853, 62881, 62909, 62938, 62965, 62993, 63031, 63058, 63086, 63114, 63142, 63171, 63198, 63226, 63254, 63282, 63310, 63338, 63366, 63394, 63422, 63450, 63478, 63506, 63534, 63562, 63590, 63618, 63646, 63674, 63712, 63738, 63766, 63794, 63822, 63851, 63878, 63906, 63934, 63962, 63990, 64018, 64046, 64074, 64103, 64130, 64158, 64186, 64214, 64242, 64270, 64309, 64338, 64366, 64395, 64422, 64450, 64478, 64506, 64534, 64562, 64590, 64618, 64647, 64674, 64703, 64731, 64759, 64787, 64815, 64853, 64882, 64910, 64938, 64966, 64994, 65022, 65050, 65078, 65107, 65134, 65162, 65191, 65218, 65247, 65275, 65303, 65331, 65369, 65398, 65426, 65454, 65482, 65510, 65538, 65566, 65594, 65623, 65650, 65679, 65707, 65735, 65763, 65802, 65830, 65858, 65898, 65925, 65953, 65994, 66022, 66050, 66078, 66106, 66134, 66162, 66190, 66219, 66246, 66274, 66302, 66330, 66358, 66386, 66414, 66442, 66471, 66498, 66526, 66555, 66583, 66611, 66639, 66666, 66695, 66722, 66751]
                xData = [-160.80, -184.05, -191.85, -198.90, -207.30, -213.60, -221.70, -227.85, -232.35, -237.75, -243.00, -248.25, -252.45, -255.90, -260.85, -265.05, -268.65, -272.70, -275.10, -279.00, -281.70, -283.20, -285.75, -287.55, -288.45, -289.95, -290.70, -291.15, -291.75, -291.15, -291.00, -289.95, -289.35, -288.00, -287.55, -286.05, -285.00, -283.05, -282.00, -279.90, -277.65, -274.50, -270.30, -266.40, -261.75, -256.05, -251.25, -247.05, -241.95, -237.30, -232.35, -227.40, -223.95, -218.55, -212.85, -207.45, -203.70, -197.85, -192.30, -186.90, -181.65, -175.65, -170.25, -163.05, -156.00, -150.60, -144.30, -138.15, -132.30, -126.45, -121.50, -116.55, -111.15, -105.15, -100.20, -95.85, -90.75, -84.60, -80.10, -75.30, -72.00, -68.10, -63.60, -59.55, -57.15, -54.90, -52.95, -50.55, -48.30, -47.40, -47.25, -47.85, -48.15, -48.90, -50.10, -53.10, -55.80, -58.95, -61.20, -64.35, -67.95, -72.15, -76.05, -80.40, -83.40, -87.75, -91.65, -94.20, -97.65, -102.00, -104.85, -109.80, -112.65, -117.90, -124.95, -131.10, -135.90, -142.35, -147.30, -153.90, -159.75, -164.70, -169.95, -174.60, -180.00, -186.75, -192.90, -197.85, -202.80, -206.40, -211.35, -214.80, -219.30, -223.20, -226.35, -229.80, -234.45, -238.20, -242.10, -244.80, -248.40, -253.20, -256.80, -260.55, -264.30, -267.30, -271.80, -274.65, -277.20, -280.50, -282.45, -285.45, -288.60, -289.50, -290.55, -291.45, -292.80, -292.35, -292.05, -292.05, -291.45, -291.45, -290.55, -287.70, -287.70, -286.80, -284.85, -281.70, -279.30, -275.85, -273.75, -269.55]
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
