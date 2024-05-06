/**
 * An extension providing a compass-bearing for a buggy located anywhere on the globe 
 * (except at the magnetic poles!), with any arbitrary mounting orientation for its microbit.
 * See the README for a detailed description of the approach, methods and algorithms.
 */

//% color=#6080e0 weight=40 icon="\uf14e" block="Heading" 
namespace heading {

    // ENUMERATIONS

    enum Dim { // Synonyms for the three magnetometer dimensions (just for brevity!)
        X = Dimension.X,
        Y = Dimension.Y,
        Z = Dimension.Z
    }

    enum View { // the axis-pairs that span the three possible magnetometer projection planes
        XY,
        YZ,
        ZX
    }

    // CONSTANTS

    const TwoPi = 2 * Math.PI
    const RadianDegrees = 360 / TwoPi
    const EnoughScanTime = 1800 // minimum acceptable scan-time
    const EnoughSamples = 70 // fewest acceptable scan samples
    const MarginalField = 30 // minimum acceptable field-strength for magnetometer readings
    const Circular = 1.03 // maximum eccentricity to consider an Ellipse as "circular"
    const LongEnough = 0.9 // for major-axis candidates, qualifying fraction of longest 

    // CLASSES

    // An Arrow is an object holding a directed vector {u,v} in both Cartesian and Polar coordinates. 
    // It also carries a time field, used to timestamp scanned samples.
    // It is used to hold a magnetometer measurement as a normalised 2D vector 
 
    class Arrow {
        u: number; // horizontal component
        v: number; // vertical component
        size: number; // polar magnitude of vector
        angle: number; // polar angle (radians anticlockwise from East)
        time: number; // (for scan samples) timestamp when this was collected

        constructor(u: number, v: number, t: number) { 
            this.u = u
            this.v = v
            this.size = Math.sqrt((u * u) + (v * v))
            this.angle = Math.atan2(v, u)
            //this.degrees = asBearing(this.angle)
            this.time = t
        }
        // cloning method
        cloneMe():Arrow {
            return new Arrow(this.u, this.v, this.time)
        }
    }

    // An Ellipse is an object holding the characteristics of the (typically) elliptical
    // view formed when projecting the scan Spin-Circle onto a 2-axis View plane.
    class Ellipse {
        plane: string; // View name (just for debug)
        uDim: number; // horizontal axis of this View
        vDim: number; // vertical axis of this View
        uOff: number; // horizontal offset needed to re-centre this Ellipse along the U-axis
        vOff: number; // vertical offset needed to re-centre this Ellipse along the V-axise circular again
    
        // calibration characteristics
        majorAxis: Arrow; // averaged major axis 
        period: number; // this View's assessment of average rotation time
        eccentricity: number; // ratio of major-axis to minor-axis magnitudes for this Ellipse
        isCircular: boolean; // flag saying this "Ellipse" View is almost circular, simplifying future handling
        fromBelow: boolean; // rotation reversal flag, reflecting this Ellipse's view of the clockwise scan



        constructor(plane: string, uDim: number, vDim: number, uOff: number, vOff: number) {
            this.plane = plane // (as a DEBUG aid)
            this.uDim = uDim
            this.vDim = vDim 
            this.uOff = uOff 
            this.vOff = vOff
            this.isCircular = false // until proved otherwise!
            this.period = -1
        }

        // Method to find the average of a set of adjacent Arrow angles (coping with cyclic angle roll-round)
        averageAngle(arrows: Arrow[]): number {
            let uSum = 0
            let vSum = 0
            // make a chain of all the arrow vectors
            for (let i = 0; i < arrows.length; i++) {
                uSum += arrows[i].u
                vSum += arrows[i].v
            }
            // the resultant vector then shows the overall direction of the chain of Arrows
            return Math.atan2(vSum, uSum)
        }

        // Method to reduce an array of seven adjacent Arrows to its central average by applying
        // the Gaussian smoothing kernel {a + 6b + 15c + 20d + 15e + 6f + g}.
        gaussianAverage(arrows: Arrow[]): Arrow {
            let uBlend = (arrows[0].u + arrows[6].u
                + 6 * (arrows[1].u + arrows[5].u)
                + 15 * (arrows[2].u + arrows[4].u)
                + 20 * arrows[3].u) / 64
            let vBlend = (arrows[0].v + arrows[6].v
                + 6 * (arrows[1].v + arrows[5].v)
                + 15 * (arrows[2].v + arrows[4].v)
                + 20 * arrows[3].v) / 64
            return new Arrow(uBlend, vBlend, arrows[3].time) // use the timestamp of the central Arrow
        }

        /// This method does two jobs:
        // 1. Finds the consensus of a set of axis-candidate Arrow angles (reversing the "opposite"
        //    ones, so they all point to the same end of the axis as the very first candidate).
        // 2. Clocks each time we pass this end, hi-jacking the "time" property of the 
        //    returned Arrow to report the average rotation period detected.
        combine(arrows: Arrow[]): Arrow {
            let turns = 0
            let endTime = 0
            let flipped = false
            let uSum = 0
            let vSum = 0
            let period = -1
            let count = arrows.length
            if (count > 0) {
                // the first candidate fixes which "end" of this axis we're choosing
                let front = arrows[0].angle
                uSum = arrows[0].u
                vSum = arrows[0].v
                let startTime = arrows[0].time
                for (let i = 1; i < count; i++) {
                    // get angle difference, as +/- 2pi
                    let deviate = ((3 * Math.PI + arrows[i].angle - front) % TwoPi) - Math.PI
                    // does it point mostly to the front? or to the back?
                    if (Math.abs(deviate) < (Math.PI/2)) {
                        // add the next arrow directly to the chain (no need to flip this one)
                        uSum += arrows[i].u
                        vSum += arrows[i].v
                        // the first unflipped Arrow after one or more flipped ones clocks a new revolution
                        if (flipped) {
                            flipped = false
                            turns++
                            endTime = arrows[i].time
                        }
                    } else { // flip this arrow before chaining it, as it's pointing the "wrong" way
                        flipped = true
                        uSum -= arrows[i].u
                        vSum -= arrows[i].v
                    }
                }
                // re-normalise the resultant's vector coordinates
                uSum /= count
                vSum /= count
                // compute the average rotation time (as long as we've made at least one complete revolution)
                if (endTime > 0) {
                    period =  (endTime - startTime) / turns
                } // else period remains at -1
            } // else period remains at -1
            // make an Arrow pointing to {uSum,vSum}, showing the overall direction of the chain of source Arrows
            return new Arrow(uSum, vSum, period)
        }
    

        // By comparing the longest and shortest radii, this method works out the eccentricity of this Ellipse.
        // It also collects possible candidates for the Ellipse major axis into this.majors[]
        // Although samples are already a rolling sum of seven readings, differentiation will amplify noise,
        // so for stability we apply an additional 7-point Gaussian smoothing as we track the first derivative.
        // We look for inflections in the slope, which occur as we pass local peaks, and push candidate values 
        // onto the majors[] axis-list for later averaging and timing analysis.
        extractAxes() {
            let arrows: Arrow[] = [] // rolling set of 7 Arrows, holding 7 adjacent samples
            let majors: Arrow[] = [] // candidate directions for major axis of Ellipse
            let longest = 0
            let shortest = 99999
            let peak = 0 // (a marker, just for debug trace)
            // collect 7 samples as Arrows and get first Gaussian sum
            for (let i = 0; i < 7; i++) {
                arrows.push(new Arrow(scanData[i][this.uDim], scanData[i][this.vDim], scanTimes[i]))
            }
            let slope: number = 99999 // marker for "first time round"
            let slopeWas: number
            let smooth = this.gaussianAverage(arrows)

            // "smooth" now contains a smoothed (w.r.t its neighbours) version of the 4th sample in this View
            let smoothWas: Arrow
            // now work through remaining samples...
            for (let i = 7; i < scanTimes.length; i++) {
                smoothWas = smooth
                arrows.shift() // drop the earliest sample & append the next one
                arrows.push(new Arrow(scanData[i][this.uDim], scanData[i][this.vDim], scanTimes[i]))
                smooth = this.gaussianAverage(arrows)
                longest = Math.max(longest, smooth.size)
                shortest = Math.min(shortest, smooth.size) 
                // now collect candidates for the major-axis angle   
                slopeWas = slope
                // differentiate to get ongoing slope
                slope = (smooth.size - smoothWas.size) / (smooth.time - smoothWas.time)
                // ensure the first two slopes always match
                if (slopeWas == 99999) slopeWas = slope
                // look for peaks where we switch from rising to falling size
                if ((slopeWas > 0) && (slope <= 0)) {
                    majors.push(smooth.cloneMe()) // copy the major axis we are passing
                    peak = 100
                }
            }

            // The ratio of the axis lengths gives the eccentricity of this Ellipse
            this.eccentricity = longest / shortest
            // Readings taken from a near-circular Ellipse won't be fore-shortened, so can skip correction!
            this.isCircular = (this.eccentricity < Circular) 

            // some collected "peaks" are only local maxima, so we must purge any whose vector length is too short
            let long = longest * LongEnough
            for (let i = 0; i < majors.length; i++) {
                if (majors[i].size < long) {
                    majors.splice(i,1)  // disqualified!
                }
            }

/* We passed each axis twice per Spin-circle revolution, so an eccentric Ellipse will produce neatly alternating
    candidates with "opposite" angles. Noisy readings mean that a more-nearly circular Ellipse may generate 
    alternating clusters of candidates as we pass each end of each axis. An almost circular Ellipse has no
    meaningful axes, but will generate multiple spurious candidates. 
    We are trying to find a good approximation to the tilt of the Ellipse's major-axis.  
    We could simply nominate the longest one detected, but while analysing rotation we average them. 
*/        
            this.majorAxis = this.combine(majors)
            this.period = this.majorAxis.time

        
            if (logging) {
                datalogger.log(
                    datalogger.createCV("view", this.plane),
                    datalogger.createCV("thetaBearing", round2(thetaBearing)),
                    datalogger.createCV("period", round2(rad2deg(this.period))),
                    datalogger.createCV("eccent.", round2(this.eccentricity)))
            }
        
        }
    
           }

    // GLOBALS

    let scanTimes: number[] = [] // sequence of time-stamps for scanned readings 
    let scanData: number[][] = [] // scanned sequence of [X,Y,Z] magnetometer readings
    let scanTime: number = 0 // duration of scan in ms
    let views: Ellipse[] = [] // the three possible elliptical views of the Spin-Circle
    let bestView = -1
    let uDim = -1 // the "horizontal" axis (called U) for the best View
    let vDim = -1 // the "vertical" axis (called V) for the best View
    let north: number // angle registered as "North", in radians counter-clockwise from U-axis
    let strength = 0 // the average magnetic field-strength observed by the magnetometer
    let fromBelow = false // set "true" if orientation means readings project backwards
    let period: number // overall assessment of average rotation time
    let rpm: number // the equivalent rotation rate in revs-per-minutetheta: number; 

    // correction parameters for future readings on bestView Ellipse
    let needsFixing: boolean // true, unless bestView Ellipse is circular
    let uOff: number // horizontal origin offset
    let vOff: number // horizontal origin offset
    let theta: number  // major-axis tilt angle (in radians anticlockwise from the U-axis)
    let thetaBearing: number; // (helpful while debugging)
    let cosTheta: number; // saved for efficiency
    let sinTheta: number; //      ditto
    let scale: number // stretch-factor for correcting foreshortened readings (= eccentricity)


    let logging = true  // logging mode flag
    let debugging = false  // test mode flag
    let dataset: string = "NONE" // test dataset to use
    let test = 0 // global selector for test-cases
    let uData: number[] = [] // U values for test cases
    let vData: number[] = [] // V values for test cases
    let northBearing = 0 // (for debug)

    // EXPORTED USER INTERFACES   

    /** 
     * Assuming the buggy is currently spinning clockwise on the spot, capture a 
     * time-stamped sequence of magnetometer readings from which to set up the compass.
     *
     * @param ms scanning-time in millisecs (long enough for more than one full rotation)    
     */

    //% block="scan for (ms) $ms" 
    //% inlineInputMode=inline 
    //% ms.shadow="timePicker" 
    //% ms.defl=0 
    //% weight=90 

    export function scan(ms: number) {
        // Every 25-30 ms over the specified duration (generally a couple of seconds),
        // magnetometer readings are sampled and a new [X,Y,Z] triple added to the scanData[] array.
        // A timestamp for each sample is also recorded in the scanTimes[] array.

        // NOTE: to smooth out jitter, each reading is always a rolling sum of SEVEN consecutive
        // readings, effectively amplifying the dynamic field range due to the Earth seven-fold 
        //(to about +/- 350 microTeslas)
        scanTimes = []
        scanData = []

        if (logging) {
            datalogger.deleteLog()
            datalogger.includeTimestamp(FlashLogTimeStampFormat.None)
        }
        if (debugging) {
            simulateScan(dataset)
            basic.pause(ms)
        } else {
            
            let index = 0
            let xRoll: number[] = []
            let yRoll: number[] = []
            let zRoll: number[] = []
            let x = 0
            let y = 0
            let z = 0
            let field = 0
            let now = input.runningTime()
            // add up the first six readings, about 25ms apart...
            for (let k = 0; k < 6; k++) {
                field = input.magneticForce(0)
                x += field
                xRoll.push(field)

                field = input.magneticForce(1)
                y += field
                yRoll.push(field)

                field = input.magneticForce(2)
                z += field
                zRoll.push(field)

                basic.pause(20)
            }

            // continue cranking out rolling sums, adding a new reading and dropping the oldest
            let finish = now + ms
            while (now < finish) {
                field = input.magneticForce(0)
                x += field
                xRoll.push(field)

                field = input.magneticForce(1)
                y += field
                yRoll.push(field)

                field = input.magneticForce(2)
                z += field
                zRoll.push(field)

                scanData.push([x, y, z])
                now = input.runningTime()
                scanTimes.push(now) 
                x -= xRoll.shift()
                y -= yRoll.shift()
                z -= zRoll.shift()
                basic.pause(20)                
            }
        }

        if (logging && !debugging) {

            datalogger.setColumnTitles("index", "t", "x", "y", "z")
            for (let i = 0; i < scanTimes.length; i++) {
                datalogger.log(
                    datalogger.createCV("index", i),
                    datalogger.createCV("t", scanTimes[i]),
                    datalogger.createCV("x", round2(scanData[i][Dim.X])),
                    datalogger.createCV("y", round2(scanData[i][Dim.Y])),
                    datalogger.createCV("z", round2(scanData[i][Dim.Z])))
            }
        }
    }



    /**
     * Analyse the scanned data to prepare for reading compass-bearings.
     * Then read the magnetometer and register the buggy's current direction as "North",
     * (i.e. the direction that will in future return zero as its heading).
     * 
     * The actual direction the buggy is pointing when this function is called could be
     * Magnetic North; or True North (compensating for local declination); or any convenient
     * direction from which to measure subsequent heading angles.
     * 
     * @returns: zero if successful, or a negative error code:
     *
     *      -1 : NOT ENOUGH SCAN DATA
     *
     *      -2 : FIELD STRENGTH TOO WEAK
     *
     *      -3 : NOT ENOUGH SCAN ROTATION
     *
     */
    //% block="set North" 
    //% inlineInputMode=inline 
    //% weight=80 
    export function setNorth(): number {
        // First analyse the scan-data to decide how best to use the magnetometer readings.
        // we'll typically need about a couple of second's worth of scanned readings...
        let nSamples = scanTimes.length
        scanTime = scanTimes[nSamples - 1] - scanTimes[0]
        if ((nSamples < EnoughSamples) || (scanTime < EnoughScanTime)) {
            return -1 // "NOT ENOUGH SCAN DATA"
        }
        // Each dimension should track a sinusoidal wave of values (generally not centred on zero).
        // The first pass finds the ranges for each axis 
        let xlo = 999
        let ylo = 999
        let zlo = 999
        let xhi = -999
        let yhi = -999
        let zhi = -999
        for (let i = 0; i < nSamples; i++) {
            xhi = Math.max(xhi, scanData[i][Dim.X])
            xlo = Math.min(xlo, scanData[i][Dim.X])
            yhi = Math.max(yhi, scanData[i][Dim.Y])
            ylo = Math.min(ylo, scanData[i][Dim.Y])
            zhi = Math.max(zhi, scanData[i][Dim.Z])
            zlo = Math.min(zlo, scanData[i][Dim.Z])
        }

        // Bail out early if the scan can't properly detect the Earth's magnetic field, perhaps due to
        // magnetic shielding, (or even because it's sitting directly over a magnetic Pole!)
        let xField = (xhi - xlo) / 2
        let yField = (yhi - ylo) / 2
        let zField = (zhi - zlo) / 2
        strength = Math.sqrt((xField * xField) + (yField * yField) + (zField * zField))
        if (strength < MarginalField) {
            return -2 // "FIELD STRENGTH TOO WEAK"
        }


        // Use the mean of these extremes as normalisation offsets
        let xOff = (xhi + xlo) / 2
        let yOff = (yhi + ylo) / 2
        let zOff = (zhi + zlo) / 2

        // apply normalisation to re-centre all of the scanData samples
        for (let i = 0; i < nSamples; i++) {
            scanData[i][Dim.X] -= xOff
            scanData[i][Dim.Y] -= yOff
            scanData[i][Dim.Z] -= zOff
        }

        // create an Ellipse instance for analysing each possible view
        views.push(new Ellipse("XY", Dim.X, Dim.Y, xOff, yOff))
        views.push(new Ellipse("YZ", Dim.Y, Dim.Z, yOff, zOff))
        views.push(new Ellipse("ZX", Dim.Z, Dim.X, zOff, xOff))

        // For each View, perform the analysis of eccentricity and Ellipse tilt-angle
        views[View.XY].extractAxes()
        views[View.YZ].extractAxes()
        views[View.ZX].extractAxes()

        // check that at least one View saw at least one complete rotation (with measurable period)...
        if (   (views[View.XY].period == -1)
            && (views[View.YZ].period == -1)
            && (views[View.ZX].period == -1) ) {
            return -3 // "NOT ENOUGH SCAN ROTATION"
        }

        // Choose the "roundest" Ellipse  --the one with lowest eccentricity.
        bestView = View.XY
        if (views[View.YZ].eccentricity < views[bestView].eccentricity) bestView = View.YZ
        if (views[View.ZX].eccentricity < views[bestView].eccentricity) bestView = View.ZX

        // periodicity is unreliable in the best View: average just the other two Views' measurements
        period = (views[0].period + views[1].period + views[2].period - views[bestView].period)/2
        rpm = 60000 / period

        // Depending on mounting orientation, the bestView Ellipse might possibly be seeing the 
        // Spin-Circle from "underneath", effectively experiencing an anti-clockwise scan.
        // Check polarity of revolution:
        fromBelow = (views[bestView].period < -1)

        // For efficiency, extract various characteristics from the bestView Ellipse
        uOff = views[bestView].uOff
        vOff = views[bestView].vOff
        scale = views[bestView].eccentricity
        theta = views[bestView].majorAxis.angle // the rotation (in radians) that aligns the major-axis clockwise with the U-axis
        thetaBearing = asBearing(theta) // (for debug)
        cosTheta = Math.cos(theta)
        sinTheta = Math.sin(theta)
        needsFixing = !views[bestView].isCircular
        // we've now finished with the scanning data, so release its memory
        scanTimes = []
        scanData = []

        // Having successfully set up the projection parameters, in the bestView get a stable fix 
        // on the current heading, which we will then designate as "North", 
        // This is the global fixed bias to be subtracted from all further readings
        north = takeSingleReading()
        northBearing = asBearing(north) // (for debug)

        if (logging) {
            datalogger.setColumnTitles("uDim", "vDim", "uOff", "vOff",
                "theta", "eccen.", "period", "north", "strength")

            datalogger.log(
                datalogger.createCV("uDim", views[bestView].uDim),
                datalogger.createCV("vDim", views[bestView].vDim),
                datalogger.createCV("uOff", round2(views[bestView].uOff)),
                datalogger.createCV("vOff", round2(views[bestView].vOff)),
                datalogger.createCV("thetaBearing", round2(thetaBearing)),
                datalogger.createCV("eccen.", round2(views[bestView].eccentricity)),
                datalogger.createCV("period", round2(views[bestView].period)),
                datalogger.createCV("north", round2(northBearing)),
                datalogger.createCV("strength", round2(strength)))
        }
        // SUCCESS!
        return 0
    }


    /**
     * Read the magnetometer
     * 
     * $returns : the current heading of the buggy (in degrees from "North")
     */
    //% block="degrees" 
    //% inlineInputMode=inline 
    //% weight=70
    export function degrees(): number {
        return asBearing(takeSingleReading() - north)
    }

    /**
     * Return average rotation time of the most recent scan 
     */
    //% block="spinTime" 
    //% inlineInputMode=inline 
    //% weight=60 
    export function scanPeriod(): number {
        if (views[bestView].period == -1) {
            return -4 // ERROR: SUCCESSFUL SCAN IS NEEDED FIRST
        } else {
            return views[bestView].period
        }
    }

    /**
     * Return average RPM of the most recent scan 
     */
    //% block="spinRPM" 
    //% inlineInputMode=inline 
    //% weight=50 
    export function spinRPM(): number {
        if (views[bestView].period <= 0) {
            return -4 // ERROR: SUCCESSFUL SCAN IS NEEDED FIRST
        } else {
            return 60000 / views[bestView].period
        }
    }

    /**
     * For scanning, wheels are rotated in opposite directions, giving a spin-rate for the 
     * selected power setting. Based on the wheel-diameter, axle-length and spin-rate, this 
     * function estimates the forward speed to be expected using that power setting.
     * (NOTE that tyre-friction when turning may make this a fairly inaccurate estimate!)
     * 
     * @param wheeelDiameter : width of wheels (in mm)
     * @param axleLength : distance betweeen mid-lines of tyres (in mm)
     * @param spinRate : the spin rotation rate reported by heading.spinRPM() for latest scan
     */
    //% block="speed using $wheelDiameter (mm), axleLength (mm) at $spinRate (RPM)" 
    //% inlineInputMode=inline 
    //% weight=50 
    export function speedUsing(wheelDiameter: number, axleLength: number, spinRate: number): number {
        if (views[bestView].period <= 0) {
            return -4 // ERROR: SUCCESSFUL SCAN IS NEEDED FIRST
        } else {
            return 99999999
        }
    }

    /**
      * Choose a test dataset to use (or "NONE")
      */
    //% block="test with: $name"
    //% inlineInputMode=inline 
    //% weight=20
    export function testDataset(name: string) {
        dataset = name
        debugging = (name.toUpperCase() != "NONE")
    }

/**
 * Choose whether to use Data Logger to grab a new test dataset into MY_DATA,
 * or to debug processing
 */
    //% block="set grab mode $on"
    //% inlineInputMode=inline 
    //% weight=10
    export function setLogMode(on: boolean) {
        logging = on
    }




    // UTILITY FUNCTIONS
    // Take the sum of seven new readings to get a stable fix on the current heading
    function takeSingleReading(): number {
        let u = 0
        let v = 0
        let reading = 0
        if (debugging) { // just choose the next test-data value
            u = uData[test]
            v = vData[test]
            test = (test + 1) % uData.length
        } else {
            // get a new sample as the sum of seven readings
            u = input.magneticForce(uDim)
            v = input.magneticForce(vDim)
            for (let i = 0; i < 6; i++) {
                basic.pause(5)
                u += input.magneticForce(uDim)
                v += input.magneticForce(vDim)
            }
        }
        // normalise the point w.r.t our Ellipse origin
        u = u - uOff
        v = v - vOff

        // Unless this Ellipse.isCircular, any {u,v} reading will be foreshortened in this View, and
        // must stretched along the Ellipse minor-axis, to place it correctly on the Spin-Circle.
        if (needsFixing) {

            // First rotate clockwise by theta, aligning the Ellipse minor-axis angle with the V-axis
            let uNew = u * cosTheta + v * sinTheta
            let vNew = v * cosTheta - u * sinTheta
            let uFixed = uNew
            // Now scale up vertically, re-balancing the axes to make the Ellipse circular
            let vFixed = vNew * scale
            // get the adjusted angle for this corrected {u,v}
            reading = Math.atan2(vFixed, uFixed)
            // finally, undo the rotation by theta
            reading += theta
        } else {
            reading = Math.atan2(v, u)
        }

        if (logging) {
            datalogger.log(
                datalogger.createCV("u", round2(u)),
                datalogger.createCV("v", round2(v)),
                datalogger.createCV("reading", Math.round(asBearing(reading))))
        }
        return reading
    }


    // convert angle from radians to degrees in the range 0..360
    function rad2deg(radians: number): number {
        let degrees = ((radians * RadianDegrees) + 360) % 360
        return degrees
    }

    // Convert angle measured in radians anticlockwise from horizontal U-axis (East)
    // to a compass bearing in degrees measured clockwise from vertical V-axis (North)
    function asBearing(angle: number): number {
        return ((90 - (angle * RadianDegrees)) + 360) % 360
    }

    // While debugging, it is useful to re-use predictable sample data for a variety of use-cases 
    // that has been captured from live runs using datalogger.
    function simulateScan(dataset: string) {
        let xData: number[] = []
        let yData: number[] = []
        let zData: number[] = []
        switch (dataset) {
            case "Xup70": // X-axis vertical; spun around X-axis; 70-degree dip
                scanTimes = [14293, 14321, 14349, 14377, 14405, 14433, 14461, 14489, 14517, 14547, 14573, 14601, 14629, 14657, 14685, 14713, 14741, 14769, 14797, 14826, 14853, 14881, 14909, 14937, 14965, 14993, 15021, 15049, 15077, 15105, 15136, 15161, 15189, 15217, 15245, 15273, 15301, 15330, 15357, 15385, 15413, 15441, 15469, 15497, 15525, 15553, 15582, 15612, 15637, 15665, 15693, 15721, 15749, 15777, 15805, 15834, 15861, 15889, 15917, 15945, 15973, 16001, 16029, 16060, 16085, 16113, 16141, 16169, 16197, 16225, 16254, 16282, 16310, 16338, 16366, 16394, 16424, 16449, 16477, 16505, 16533, 16561, 16590, 16618, 16646, 16674, 16702, 16730, 16761, 16789, 16817, 16845, 16873, 16901, 16930, 16958, 16986, 17014, 17042, 17070, 17101, 17129, 17157, 17185, 17213, 17241, 17270, 17301, 17329, 17357, 17385, 17413, 17441, 17470, 17498, 17526, 17554, 17585, 17613, 17641, 17669, 17697, 17726, 17754, 17782, 17813, 17841, 17870, 17897, 17925, 17954, 17982, 18010, 18041, 18069, 18097, 18126, 18154, 18182, 18213, 18241, 18269, 18298, 18326, 18354, 18385, 18413, 18442, 18470, 18498, 18526, 18557, 18585, 18614, 18642, 18670, 18701, 18729, 18758, 18786, 18817, 18845, 18874, 18902, 18930, 18961, 18990, 19018, 19046, 19077, 19106, 19134]
                xData = [146.85, 147.30, 146.25, 146.25, 146.85, 147.30, 147.00, 147.90, 147.45, 148.50, 148.50, 149.70, 149.55, 150.30, 150.15, 150.30, 150.30, 152.10, 151.65, 153.15, 153.00, 153.00, 154.20, 154.05, 153.45, 152.70, 152.10, 152.85, 153.30, 152.70, 153.30, 153.30, 154.20, 153.90, 152.55, 152.40, 153.15, 153.15, 153.45, 152.85, 153.00, 153.15, 152.70, 151.65, 150.90, 149.70, 149.40, 148.20, 148.50, 148.20, 149.70, 150.45, 151.05, 151.05, 151.20, 150.90, 149.85, 147.75, 146.25, 145.80, 145.80, 145.65, 145.50, 145.95, 145.35, 145.35, 144.45, 144.30, 145.20, 145.95, 146.10, 145.95, 145.50, 146.25, 145.80, 145.05, 143.70, 143.25, 143.70, 143.85, 143.25, 142.65, 143.10, 142.95, 142.35, 142.80, 142.65, 143.55, 144.45, 144.75, 145.50, 147.15, 146.70, 147.45, 147.00, 146.85, 146.10, 145.35, 145.65, 146.10, 145.95, 145.65, 146.40, 147.45, 148.35, 147.45, 147.75, 146.85, 147.90, 147.60, 147.00, 147.00, 146.70, 147.15, 148.50, 147.90, 147.45, 147.45, 147.60, 149.70, 148.65, 149.10, 149.25, 150.00, 150.90, 151.80, 151.95, 152.40, 152.25, 151.95, 151.95, 151.05, 151.05, 149.40, 150.60, 150.75, 152.55, 152.40, 152.70, 152.55, 153.75, 153.45, 154.65, 154.50, 154.95, 155.40, 154.80, 154.20, 154.35, 153.30, 151.95, 151.50, 151.50, 151.50, 151.65, 151.80, 150.90, 151.95, 151.95, 152.10, 151.95, 151.05, 150.45, 150.45, 149.40, 149.40, 148.65, 148.35, 149.10, 148.80, 148.05, 148.35]
                yData = [219.75, 213.60, 206.70, 199.80, 194.70, 189.15, 184.80, 181.95, 178.05, 175.65, 172.80, 170.70, 169.35, 167.70, 166.20, 166.50, 168.15, 170.40, 171.90, 173.70, 175.50, 177.45, 179.55, 181.95, 184.95, 187.95, 191.70, 196.35, 201.00, 206.40, 211.20, 216.30, 222.30, 228.45, 234.30, 240.75, 247.20, 253.65, 260.70, 266.85, 273.00, 280.65, 288.45, 294.75, 301.80, 310.35, 318.00, 325.95, 333.30, 340.95, 349.95, 358.80, 366.45, 374.85, 383.40, 392.40, 400.20, 406.80, 412.80, 417.30, 421.50, 424.35, 425.40, 425.70, 426.30, 425.10, 424.05, 422.85, 422.25, 420.75, 419.25, 417.90, 415.95, 413.55, 410.85, 408.15, 405.45, 402.00, 397.05, 393.30, 388.50, 383.70, 377.85, 372.30, 366.90, 361.65, 354.90, 349.95, 343.50, 336.90, 329.85, 323.10, 316.35, 311.40, 304.65, 298.65, 292.65, 287.25, 281.85, 275.55, 268.50, 261.00, 254.85, 248.70, 241.80, 235.20, 229.35, 224.25, 219.75, 214.50, 210.15, 205.80, 201.45, 198.45, 193.95, 189.90, 186.30, 182.25, 178.50, 175.35, 171.15, 168.75, 167.55, 166.05, 164.85, 163.50, 163.80, 164.55, 165.30, 165.30, 167.40, 169.65, 172.80, 174.45, 176.70, 180.45, 184.35, 187.65, 191.70, 196.35, 200.40, 206.40, 210.75, 215.85, 220.65, 226.65, 232.80, 240.00, 246.45, 253.80, 261.75, 270.75, 280.20, 288.00, 295.95, 303.60, 311.70, 319.50, 327.00, 332.70, 339.45, 346.05, 353.25, 359.10, 365.70, 371.70, 378.30, 383.85, 389.70, 394.20, 399.45, 404.25]
                zData = [-179.40, -174.45, -168.75, -162.75, -155.70, -148.50, -141.30, -134.70, -126.00, -117.75, -109.50, -101.40, -93.75, -86.10, -78.90, -72.15, -64.65, -58.20, -51.60, -45.30, -39.30, -34.05, -27.90, -23.55, -18.90, -13.95, -8.55, -3.00, 3.00, 7.35, 12.00, 16.65, 21.15, 24.00, 26.70, 28.20, 29.85, 31.35, 32.55, 32.25, 33.00, 32.55, 33.00, 31.95, 31.65, 30.45, 28.95, 25.95, 23.70, 19.95, 16.50, 10.50, 3.60, -2.25, -10.50, -19.35, -28.80, -39.45, -49.20, -58.35, -69.15, -76.80, -85.80, -93.90, -100.95, -106.50, -112.65, -116.70, -121.95, -126.30, -131.25, -135.00, -141.90, -147.30, -153.45, -157.50, -162.45, -166.20, -172.35, -176.25, -180.45, -185.25, -190.95, -195.00, -197.85, -199.95, -202.50, -204.75, -206.10, -207.45, -207.60, -210.00, -210.75, -211.95, -211.20, -210.60, -209.85, -209.10, -207.15, -205.35, -203.70, -202.50, -201.15, -197.55, -196.35, -192.15, -188.25, -183.00, -178.35, -173.55, -169.80, -164.25, -161.25, -157.05, -152.40, -147.90, -142.20, -136.95, -131.10, -123.75, -116.10, -109.50, -102.00, -94.80, -87.00, -80.25, -75.45, -70.35, -64.35, -58.20, -53.70, -48.90, -43.20, -35.70, -29.40, -23.55, -18.30, -10.95, -5.85, -1.20, 2.10, 6.75, 11.10, 14.55, 17.10, 21.30, 25.20, 27.75, 29.25, 30.15, 32.10, 32.85, 32.55, 31.95, 33.00, 31.95, 30.45, 28.05, 26.40, 24.75, 22.35, 19.20, 16.50, 13.80, 9.30, 4.50, -0.90, -7.05, -13.95, -21.30, -28.95, -36.75]
                uData = [417.00, 414.15, 339.15, 238.80, 172.35, 171.90, 248.40, 348.00, 413.25]
                vData = [-53.10, -144.45, -204.15, -194.40, -115.05, -29.55, 31.80, 20.55, -52.35]
                break

            case "Yup70": // Y-axis vertical; spun around Y-axis; 70-degree dip
                scanTimes = [10844, 10872, 10900, 10928, 10956, 10984, 11012, 11040, 11068, 11098, 11124, 11152, 11180, 11208, 11236, 11264, 11292, 11320, 11348, 11376, 11404, 11432, 11460, 11488, 11516, 11544, 11572, 11600, 11628, 11659, 11684, 11712, 11740, 11768, 11796, 11824, 11852, 11880, 11908, 11936, 11964, 11992, 12020, 12048, 12076, 12107, 12132, 12160, 12188, 12216, 12244, 12272, 12300, 12328, 12356, 12384, 12413, 12440, 12468, 12496, 12524, 12552, 12583, 12608, 12636, 12665, 12692, 12720, 12748, 12776, 12804, 12832, 12860, 12888, 12917, 12947, 12972, 13000, 13028, 13056, 13084, 13112, 13140, 13168, 13196, 13224, 13252, 13280, 13311, 13336, 13364, 13392, 13421, 13448, 13476, 13504, 13532, 13560, 13588, 13616, 13647, 13672, 13700, 13728, 13756, 13784, 13812, 13844, 13872, 13900, 13928, 13956, 13984, 14012, 14041, 14069, 14097, 14127, 14152, 14180, 14208, 14236, 14264, 14293, 14320, 14352, 14380, 14408, 14436, 14464, 14492, 14520, 14552, 14580, 14608, 14636, 14665, 14692, 14721, 14752, 14780, 14808, 14837, 14865, 14893, 14924, 14953, 14980, 15008, 15036, 15068, 15096, 15124, 15152, 15180, 15212, 15240, 15268, 15296, 15325, 15356, 15384, 15412, 15441, 15472, 15500, 15529, 15557, 15588, 15616, 15644, 15676]
                xData = [-223.35, -253.35, -259.80, -264.90, -269.40, -274.35, -278.25, -280.35, -282.45, -283.50, -284.40, -285.45, -286.65, -287.70, -288.15, -289.20, -289.65, -289.65, -289.05, -288.15, -286.35, -285.00, -283.80, -281.70, -279.30, -277.95, -274.05, -270.60, -267.60, -263.55, -259.65, -256.20, -252.30, -249.30, -245.25, -241.05, -236.25, -231.60, -226.80, -221.85, -217.65, -213.90, -209.40, -205.20, -201.15, -196.80, -191.55, -186.30, -180.45, -175.05, -169.95, -163.65, -158.70, -152.70, -147.30, -142.05, -137.70, -132.15, -127.20, -121.80, -116.70, -110.85, -105.15, -98.40, -93.60, -88.35, -82.05, -77.40, -73.95, -70.05, -66.45, -62.10, -59.10, -57.00, -54.75, -51.90, -50.25, -48.45, -48.00, -47.10, -46.35, -45.90, -45.15, -44.85, -46.05, -45.30, -47.10, -48.75, -50.55, -51.75, -52.95, -54.90, -58.05, -60.30, -63.15, -66.15, -69.90, -75.00, -79.35, -84.30, -89.70, -94.20, -98.85, -104.55, -109.80, -115.80, -120.90, -126.00, -132.15, -139.35, -145.65, -151.65, -156.30, -161.25, -166.35, -170.85, -175.20, -180.30, -184.95, -190.05, -195.45, -200.70, -205.80, -210.90, -214.50, -219.15, -223.05, -227.70, -231.60, -236.25, -240.45, -245.40, -250.50, -255.00, -259.05, -262.95, -266.55, -269.55, -272.25, -273.30, -275.40, -277.65, -280.35, -282.75, -284.55, -286.05, -288.15, -289.35, -290.85, -291.45, -290.55, -290.85, -290.40, -289.80, -290.10, -287.85, -285.75, -285.00, -283.35, -282.30, -280.80, -277.65, -276.00, -273.60, -271.20, -267.75, -263.10, -259.50, -254.70, -249.60, -244.35, -237.60]
                yData = [-27.30, -27.00, -26.85, -26.70, -26.40, -26.85, -27.00, -26.55, -26.55, -27.00, -27.30, -28.05, -27.60, -27.60, -27.30, -27.00, -26.40, -27.15, -26.85, -27.45, -27.90, -28.65, -29.10, -29.70, -29.40, -30.00, -29.85, -30.30, -29.25, -28.95, -28.50, -28.20, -27.30, -27.15, -26.85, -27.75, -28.35, -28.50, -28.35, -28.20, -28.20, -28.95, -29.10, -29.40, -30.30, -30.60, -31.35, -31.95, -31.20, -31.05, -30.45, -30.00, -30.00, -29.25, -29.25, -29.25, -28.95, -29.40, -28.80, -28.80, -29.85, -29.85, -29.25, -30.00, -30.30, -30.30, -29.70, -30.00, -30.90, -32.10, -32.40, -31.95, -33.00, -34.50, -34.35, -33.00, -33.30, -33.30, -33.30, -33.00, -32.70, -32.55, -33.30, -32.70, -32.55, -33.15, -33.45, -34.20, -33.90, -34.35, -34.20, -33.45, -33.30, -33.00, -31.65, -31.95, -30.90, -30.30, -30.00, -29.40, -29.85, -29.85, -28.95, -28.50, -28.20, -28.20, -29.10, -28.50, -28.05, -28.65, -29.10, -29.40, -29.85, -28.95, -28.95, -29.10, -28.50, -27.75, -28.05, -27.00, -27.30, -26.10, -25.35, -25.05, -25.20, -24.90, -24.75, -23.55, -24.15, -24.60, -24.75, -25.20, -24.90, -26.40, -27.00, -27.45, -27.45, -28.35, -28.35, -28.80, -27.75, -28.20, -28.20, -29.25, -27.90, -27.75, -27.30, -28.05, -27.30, -27.60, -27.75, -28.20, -28.05, -27.90, -27.60, -28.35, -27.30, -26.70, -27.15, -27.90, -28.50, -28.05, -27.75, -27.90, -27.30, -26.70, -26.55, -26.55, -27.15, -27.60, -28.05, -28.50]
                zData = [-131.85, -148.35, -143.10, -137.10, -131.70, -127.80, -121.20, -115.80, -109.50, -103.05, -96.45, -90.00, -83.40, -79.05, -73.20, -67.65, -62.40, -57.45, -53.55, -47.25, -41.10, -34.80, -30.60, -25.80, -20.25, -13.50, -9.30, -4.35, 0.30, 6.30, 11.10, 14.70, 18.15, 22.95, 26.55, 29.40, 31.35, 34.50, 37.20, 40.05, 41.55, 42.90, 45.30, 47.40, 48.75, 50.25, 52.20, 52.95, 54.30, 52.95, 52.80, 53.25, 53.25, 51.90, 51.45, 50.40, 50.85, 49.95, 48.15, 46.50, 43.95, 41.85, 38.10, 34.05, 29.70, 24.75, 20.70, 16.65, 12.15, 8.25, 3.75, -1.65, -6.15, -12.75, -18.90, -25.65, -32.85, -38.40, -43.35, -49.05, -54.00, -59.40, -64.35, -68.40, -73.95, -79.35, -84.45, -89.40, -93.90, -99.45, -105.15, -110.25, -114.75, -120.00, -126.15, -130.05, -135.60, -140.10, -144.60, -148.50, -153.00, -156.90, -162.00, -164.55, -166.80, -169.95, -174.15, -177.00, -178.80, -180.75, -182.85, -185.40, -186.60, -186.45, -185.85, -187.20, -187.20, -186.45, -185.55, -185.40, -185.25, -184.05, -181.05, -179.85, -178.50, -177.00, -174.30, -172.20, -171.00, -168.45, -165.15, -161.25, -157.50, -153.60, -149.25, -144.75, -141.00, -138.60, -135.75, -130.95, -127.35, -123.90, -119.55, -114.90, -108.75, -102.90, -98.40, -91.95, -85.80, -80.40, -75.15, -69.75, -63.30, -57.45, -52.80, -47.85, -42.60, -36.90, -31.50, -28.80, -23.55, -18.90, -13.80, -8.85, -4.35, 0.60, 6.45, 10.20, 15.15, 19.80, 25.20, 28.20]
                uData = [-22.65, -111.6, -175.5, -180.15, -115.35, -20.55, 44.70, 43.95, -17.10]
                vData = [-54.75, -57.3, -127.95, -219.3, -286.8, -280.20, -210.90, -119.40, -52.95]
                break

            case "Zdn70": // Z-axis vertically down; spun around Z-axis; 70-degree dip
                scanTimes = [41025, 41053, 41081, 41109, 41137, 41165, 41193, 41221, 41249, 41279, 41305, 41333, 41361, 41389, 41417, 41445, 41473, 41501, 41529, 41558, 41585, 41613, 41641, 41669, 41697, 41725, 41753, 41781, 41810, 41840, 41865, 41893, 41921, 41949, 41977, 42005, 42033, 42062, 42089, 42117, 42145, 42173, 42201, 42229, 42257, 42285, 42316, 42341, 42369, 42397, 42425, 42453, 42481, 42509, 42537, 42566, 42593, 42621, 42650, 42678, 42706, 42734, 42764, 42789, 42818, 42845, 42873, 42901, 42929, 42957, 42985, 43013, 43041, 43070, 43098, 43129, 43157, 43185, 43213, 43241, 43269, 43297, 43325, 43353, 43382, 43410, 43438, 43466, 43497, 43525, 43553, 43581, 43609, 43637, 43665, 43694, 43722, 43750, 43778, 43809, 43837, 43866, 43893, 43921, 43949, 43978, 44005, 44037, 44065, 44093, 44121, 44149, 44178, 44206, 44234, 44262, 44293, 44321, 44350, 44377, 44406, 44434, 44462, 44490, 44518, 44549, 44577, 44605, 44634, 44662, 44690, 44718, 44749, 44777, 44805, 44833, 44862, 44890, 44918, 44949, 44977, 45005, 45034, 45062, 45090, 45121, 45149, 45178, 45206, 45234, 45265, 45293, 45322, 45350, 45378, 45409, 45437, 45466, 45494, 45522, 45553, 45581, 45610, 45638, 45669, 45697, 45726, 45754, 45785, 45813, 45842, 45870]
                xData = [-286.35, -285.45, -284.55, -283.35, -280.35, -276.15, -272.25, -268.05, -264.3, -259.2, -254.1, -249.3, -244.65, -238.65, -233.25, -227.25, -219.6, -211.8, -202.5, -193.2, -184.5, -175.05, -164.85, -155.4, -145.8, -136.95, -128.85, -119.7, -110.1, -101.85, -95.1, -88.05, -81.45, -74.25, -68.4, -64.2, -58.5, -53.55, -50.4, -48.45, -45.9, -43.65, -42.3, -42, -42, -41.25, -40.05, -41.1, -41.25, -42.15, -43.35, -44.4, -45.45, -46.95, -49.05, -52.05, -55.35, -59.25, -63, -68.4, -75.15, -80.4, -87.15, -92.85, -98.7, -105.45, -110.7, -115.8, -122.4, -127.2, -133.35, -139.35, -145.05, -151.5, -157.05, -162.45, -169.05, -174.75, -180.45, -186.75, -193.2, -200.7, -206.7, -212.55, -218.4, -224.4, -229.2, -234.45, -238.2, -243.15, -247.5, -251.25, -254.1, -258, -260.85, -264.3, -267.6, -270.3, -273.45, -276.15, -279.15, -281.7, -282.9, -283.8, -285.6, -286.65, -287.7, -286.8, -285.9, -285, -283.5, -281.85, -279.9, -277.05, -274.65, -272.4, -269.7, -266.7, -262.95, -258.75, -254.4, -250.65, -245.25, -240.75, -235.95, -230.4, -224.55, -219.15, -213, -206.7, -200.4, -193.05, -185.85, -178.8, -172.2, -164.1, -157.05, -148.35, -141.3, -133.65, -126, -117.45, -109.8, -101.85, -94.95, -87.75, -82.05, -75.15, -69.75, -64.5, -60, -55.65, -52.35, -48.6, -47.1, -44.7, -43.8, -42.9, -42.75, -42.6, -44.25, -45.75, -48.6, -51.6, -54.9, -60, -65.55, -70.5, -75.9, -82.05, -88.65, -96.3]
                yData = [324.75, 334.05, 343.05, 352.2, 360.9, 369.15, 375.9, 383.7, 390, 396.15, 400.8, 406.05, 411.75, 417.15, 421.05, 426.15, 430.2, 434.25, 437.25, 438.45, 440.7, 442.95, 442.2, 441.45, 440.55, 437.7, 435.75, 432.75, 428.55, 425.7, 421.35, 416.4, 411.15, 405.45, 398.4, 391.35, 382.95, 374.85, 366.6, 359.55, 351.75, 344.4, 337.05, 329.85, 323.7, 316.8, 310.5, 304.5, 298.05, 291.75, 285.9, 279.3, 272.4, 266.55, 259.2, 252.6, 245.85, 238.8, 233.25, 227.7, 221.25, 217.05, 212.7, 208.8, 205.5, 201.45, 199.2, 196.95, 195, 193.35, 192.15, 191.4, 190.35, 189.9, 189.9, 188.55, 189, 190.35, 191.25, 193.35, 194.25, 195.75, 198.75, 202.2, 204.15, 206.4, 209.55, 215.55, 219.15, 223.35, 226.2, 230.1, 235.35, 238.95, 241.2, 245.55, 250.8, 256.2, 261, 266.4, 272.7, 279, 286.5, 293.4, 300.3, 308.85, 316.65, 324.75, 332.4, 339.3, 346.65, 354.3, 360.15, 366, 371.85, 378.15, 385.2, 390.9, 396.3, 402.15, 406.95, 411.9, 417.15, 420.3, 422.85, 426, 429, 432.45, 434.7, 436.95, 439.2, 442.2, 442.65, 444.15, 444, 444.75, 444.15, 443.55, 441.45, 440.7, 437.55, 435.3, 430.8, 426.6, 421.95, 415.95, 410.1, 403.8, 396.75, 389.25, 381, 372.75, 365.55, 356.25, 346.65, 337.2, 328.5, 319.2, 308.85, 297.75, 287.55, 278.85, 269.85, 261, 251.1, 243, 235.95, 229.05, 221.4, 215.25, 209.7, 206.4]
                zData = [209.4, 209.1, 208.95, 209.25, 207.75, 206.4, 205.2, 204.6, 204, 203.1, 202.05, 202.95, 202.65, 202.5, 202.5, 202.05, 202.5, 203.25, 203.1, 203.1, 201.9, 201.75, 200.55, 199.5, 198.9, 199.35, 198.75, 198.9, 199.65, 201.6, 202.05, 201.75, 201.3, 202.95, 204.45, 204, 204.6, 206.25, 207.6, 209.25, 209.85, 210.45, 211.5, 212.1, 212.25, 213.6, 213.45, 214.5, 215.4, 215.1, 215.55, 216.6, 215.85, 216.6, 216.3, 217.35, 219, 219.45, 218.4, 219.45, 219.75, 220.2, 218.7, 219.15, 219.15, 220.95, 221.1, 221.55, 223.35, 225.3, 224.85, 225, 225.15, 225.75, 226.5, 224.55, 223.05, 224.25, 224.1, 223.8, 222.9, 222.45, 222, 222.3, 220.35, 221.1, 220.35, 220.35, 219.6, 219.9, 219.45, 220.5, 217.8, 217.65, 217.2, 217.05, 216.45, 216, 214.65, 215.1, 214.35, 212.85, 211.65, 212.25, 211.2, 209.7, 208.95, 208.2, 207.75, 206.55, 204, 203.4, 203.7, 204, 202.5, 202.35, 201.75, 202.35, 201.75, 201.15, 201, 202.2, 201.6, 201.45, 200.55, 201.45, 200.7, 200.25, 198.75, 200.1, 200.1, 200.55, 199.2, 200.85, 200.55, 201.6, 199.95, 199.8, 200.25, 201.6, 201, 201.15, 200.85, 201.75, 202.65, 201.9, 202.65, 203.4, 203.85, 204.3, 205.95, 207.9, 210.45, 210.45, 210.9, 212.4, 214.05, 214.65, 213.9, 214.35, 215.1, 215.7, 216.45, 217.05, 217.35, 219.3, 218.85, 219.45, 220.8, 221.1, 221.55, 222]
                uData = [-47.85, -49.35, -118.35, -211.95, -280.95, -275.85, -208.2, -116.55, -46.05]
                vData = [357.6, 261.9, 197.25, 205.05, 273.75, 375.6, 436.35, 433.05, 361.95]
                break

            case "strange": // Z-axis vertically down; spun around Z-axis; 70-degree dip
                scanTimes = [31809, 31837, 31865, 31893, 31921, 31949, 31977, 32005, 32033, 32063, 32089, 32117, 32145, 32173, 32201, 32229, 32257, 32286, 32313, 32341, 32369, 32397, 32425, 32453, 32481, 32509, 32537, 32565, 32593, 32624, 32649, 32677, 32705, 32734, 32761, 32789, 32817, 32845, 32873, 32901, 32929, 32957, 32985, 33013, 33041, 33069, 33100, 33125, 33153, 33181, 33209, 33237, 33265, 33293, 33321, 33349, 33377, 33405, 33434, 33461, 33490, 33517, 33548, 33573, 33601, 33629, 33657, 33685, 33713, 33741, 33769, 33798, 33826, 33854, 33882, 33912, 33937, 33965, 33993, 34021, 34049, 34078, 34106, 34134, 34162, 34190, 34218, 34249, 34277, 34305, 34333, 34361, 34389, 34418, 34446, 34474, 34502, 34530, 34558, 34589, 34617, 34645, 34673, 34701, 34729, 34758, 34786, 34817, 34845, 34874, 34901, 34929, 34958, 34986, 35014, 35042, 35070, 35101, 35129, 35157, 35186, 35214, 35242, 35270, 35298, 35329, 35357, 35385, 35414, 35442, 35470, 35498, 35529, 35557, 35585, 35614, 35642, 35670, 35698, 35729, 35757, 35786, 35814, 35842, 35870, 35901, 35929, 35958, 35986, 36014, 36045, 36073, 36102, 36130, 36158, 36189, 36217, 36246, 36274, 36305, 36334, 36362, 36390, 36418, 36449, 36478, 36506, 36537, 36565, 36594, 36622, 36653]
                xData = [-421.35, -428.85, -433.95, -441, -447, -454.2, -460.95, -466.65, -470.55, -476.4, -481.8, -486.3, -491.55, -495.9, -499.5, -504.15, -507.75, -510.9, -513.9, -514.95, -517.05, -517.35, -518.25, -518.55, -518.55, -519.9, -519.9, -520.05, -521.7, -521.85, -521.55, -520.35, -517.35, -515.25, -512.4, -509.25, -505.5, -501.6, -498, -494.25, -490.8, -486.9, -481.65, -476.7, -472.5, -467.55, -462.6, -457.65, -451.05, -446.4, -441.15, -433.8, -428.4, -421.5, -414.75, -407.4, -400.2, -392.7, -387, -378.6, -372.9, -366.15, -360.15, -353.55, -347.25, -340.95, -336.15, -330.9, -325.35, -319.95, -315, -310.5, -306.15, -301.95, -297.75, -294.3, -292.8, -290.7, -288, -286.5, -285.45, -284.25, -283.95, -283.2, -283.2, -284.25, -285.9, -286.65, -289.5, -292.65, -296.7, -300.45, -306, -309.75, -316.05, -321.6, -326.7, -332.85, -339.9, -346.05, -354.45, -360.9, -367.95, -375.3, -382.05, -388.95, -395.1, -400.8, -408.15, -413.85, -419.85, -425.55, -431.1, -437.25, -442.8, -448.5, -453.75, -459.3, -464.1, -469.35, -474.15, -478.35, -481.65, -485.25, -488.55, -492.75, -495.45, -498.15, -501.75, -504.45, -508.05, -510.6, -513.3, -515.25, -517.65, -518.7, -520.95, -522.6, -523.5, -523.95, -523.65, -523.5, -524.55, -522.75, -521.1, -519.9, -517.65, -516.3, -513.45, -509.25, -506.1, -502.05, -497.25, -493.35, -489, -483.6, -477.3, -471.9, -465, -457.5, -448.95, -440.25, -431.85, -423.9, -414.15, -405.75, -396.3, -387.6, -377.7, -369.9, -361.2, -353.85]
                yData = [202.05, 202.35, 203.4, 204.9, 205.5, 206.85, 208.8, 210.45, 212.1, 213.45, 214.8, 217.8, 221.4, 225, 228.75, 232.95, 238.35, 243.15, 248.1, 251.85, 256.35, 261.9, 267.75, 272.85, 279.15, 284.7, 290.55, 295.05, 299.7, 304.5, 310.8, 316.5, 323.1, 328.65, 336, 343.2, 349.95, 355.35, 360.9, 365.55, 371.1, 375.45, 379.8, 384.15, 388.5, 393.15, 397.5, 401.1, 404.85, 407.85, 411.15, 414.75, 416.4, 418.2, 419.85, 422.7, 423.9, 425.1, 424.8, 425.25, 425.55, 426.3, 424.95, 424.2, 421.35, 419.25, 417.15, 413.85, 410.7, 407.7, 404.25, 401.55, 397.5, 393.9, 390.3, 385.95, 381.15, 376.5, 370.65, 365.4, 358.2, 350.85, 343.65, 335.7, 328.05, 320.25, 311.7, 303.9, 296.25, 288.45, 281.4, 274.35, 268.35, 262.05, 256.65, 251.1, 245.55, 239.55, 233.55, 228.45, 224.7, 220.2, 216.9, 213.45, 211.5, 209.25, 206.4, 204.45, 203.55, 202.35, 201.9, 201.45, 201, 202.5, 203.55, 204, 204.75, 206.4, 207.6, 209.85, 210.9, 212.4, 214.8, 218.25, 221.1, 223.95, 227.7, 231.75, 234.9, 238.65, 242.7, 246.9, 251.85, 256.05, 261.3, 267, 271.8, 276.15, 279.9, 284.1, 289.65, 294.15, 299.7, 305.4, 310.65, 318, 325.2, 331.2, 337.5, 343.2, 349.35, 357.15, 362.25, 367.2, 373.5, 379.65, 385.2, 390.45, 395.25, 401.1, 405.6, 409.05, 412.2, 415.35, 418.8, 420.15, 421.2, 423, 423.9, 424.05, 423.6, 422.7]
                zData = [454.5, 451.2, 448.2, 444.15, 440.55, 436.8, 432.6, 428.85, 426, 422.25, 418.8, 415.05, 411.15, 407.85, 404.4, 400.65, 395.4, 391.35, 387.15, 383.55, 379.2, 375, 370.65, 367.35, 363, 359.1, 354.9, 350.55, 348, 344.4, 341.85, 338.85, 336.45, 334.05, 331.95, 328.65, 327.45, 325.65, 324.45, 323.7, 321.9, 321.6, 321, 319.65, 320.1, 319.8, 320.1, 321.15, 321.75, 322.95, 324, 325.05, 327.6, 329.1, 330.45, 333, 335.85, 339.3, 341.7, 344.25, 348, 351.6, 355.5, 357.9, 361.8, 365.85, 369.75, 373.2, 376.95, 381, 385.2, 389.1, 392.55, 396.75, 400.8, 405.75, 409.65, 413.7, 418.8, 424.2, 428.25, 432.3, 436.05, 440.4, 445.35, 448.5, 452.1, 455.85, 459, 462.9, 465.6, 467.7, 469.05, 470.85, 472.35, 474.15, 473.7, 473.25, 473.4, 473.55, 472.2, 471.3, 469.2, 468.15, 466.5, 463.95, 462.3, 461.55, 459.75, 457.65, 455.4, 452.55, 450.6, 447.9, 444.75, 441.9, 440.25, 437.25, 434.1, 430.5, 426.75, 423.45, 419.4, 414.6, 411, 408.45, 405.75, 402.45, 397.2, 393.75, 390.15, 386.25, 381.45, 376.65, 372.3, 369.75, 366.15, 363, 360.3, 357.75, 355.05, 352.8, 350.1, 346.65, 343.95, 340.5, 337.95, 334.2, 331.95, 328.65, 326.7, 324.45, 322.8, 322.2, 322.95, 321.6, 322.05, 321.9, 322.05, 323.7, 323.4, 323.1, 325.5, 327, 329.85, 332.25, 334.95, 338.7, 343.95, 347.55, 351.75, 356.25]
                uData = [-280.35, -296.1, -372.75, -466.8, -520.2, -503.7, -424.65, -336.45, -281.7]
                vData = [358.65, 277.65, 214.05, 207.15, 267.15, 352.35, 415.95, 420, 361.35]
                break

            case "distort": // Z-axis down, but local magnetic distortion on buggy
                scanTimes = [63117, 63145, 63173, 63201, 63229, 63257, 63285, 63313, 63341, 63371, 63397, 63425, 63453, 63481, 63509, 63537, 63565, 63593, 63621, 63649, 63677, 63705, 63733, 63761, 63789, 63817, 63845, 63873, 63901, 63932, 63957, 63985, 64013, 64041, 64069, 64097, 64125, 64153, 64181, 64209, 64237, 64265, 64294, 64321, 64349, 64377, 64408, 64433, 64461, 64489, 64517, 64545, 64573, 64601, 64629, 64657, 64685, 64713, 64742, 64770, 64798, 64826, 64856, 64881, 64909, 64937, 64965, 64994, 65021, 65049, 65077, 65105, 65133, 65161, 65189, 65220, 65245, 65273, 65301, 65329, 65357, 65386, 65414, 65442, 65470, 65498, 65526, 65554, 65584, 65609, 65637, 65665, 65693, 65721, 65749, 65778, 65806, 65834, 65862, 65893, 65921, 65949, 65977, 66006, 66033, 66062, 66090, 66121, 66149, 66177, 66205, 66233, 66262, 66290, 66318, 66346, 66377, 66405, 66433, 66462, 66490, 66518, 66546, 66574, 66605, 66633, 66661, 66690, 66718, 66746, 66774, 66802, 66833, 66861, 66889, 66918, 66946, 66974, 67005, 67033, 67061, 67090, 67118, 67146, 67177, 67205, 67233, 67261, 67290, 67318, 67349, 67377, 67406, 67434, 67462, 67493, 67522, 67549, 67578, 67606, 67637, 67665, 67694, 67722, 67753, 67781, 67810, 67838, 67869, 67898, 67926, 67957]
                xData = [824.7, 822.6, 821.7, 820.65, 820.65, 820.65, 822.3, 825.45, 828.6, 832.35, 836.85, 842.7, 848.4, 854.1, 860.7, 867.75, 875.7, 883.65, 890.4, 899.25, 908.55, 918.6, 928.5, 939, 947.7, 958.65, 970.05, 980.4, 989.85, 1000.2, 1009.35, 1020.6, 1030.8, 1040.1, 1048.35, 1058.7, 1068, 1077, 1087.8, 1097.7, 1107.3, 1118.25, 1126.5, 1135.5, 1144.65, 1152.3, 1159.5, 1167.45, 1175.7, 1182.9, 1188.75, 1195.2, 1201.2, 1207.5, 1211.55, 1215.6, 1220.25, 1225.8, 1230.15, 1234.35, 1238.55, 1243.65, 1247.7, 1252.2, 1255.05, 1258.05, 1260.45, 1261.95, 1263.6, 1263, 1263.6, 1263.6, 1262.7, 1260.9, 1259.55, 1256.55, 1255.35, 1251, 1248, 1244.4, 1241.1, 1236.75, 1231.5, 1226.55, 1221.9, 1216.2, 1211.4, 1205.85, 1200.75, 1196.1, 1190.1, 1185, 1179.6, 1173.3, 1166.85, 1161.15, 1154.85, 1149.75, 1143, 1137, 1129.5, 1123.65, 1115.7, 1108.2, 1100.1, 1092.45, 1084.65, 1077.3, 1068.9, 1061.7, 1055.25, 1048.65, 1042.65, 1036.65, 1031.85, 1027.35, 1021.95, 1016.1, 1010.25, 1005.3, 999.75, 994.5, 988.5, 983.4, 977.7, 971.55, 965.25, 958.35, 949.65, 941.4, 932.25, 924, 916.5, 907.8, 901.35, 895.05, 888.3, 883.8, 879, 873.3, 868.65, 863.7, 860.1, 856.8, 852, 847.8, 845.55, 842.1, 837.3, 834.45, 831.75, 830.85, 829.95, 828, 828, 829.65, 829.65, 831.45, 832.5, 835.2, 838.95, 842.25, 846.6, 851.55, 858.15, 864.75, 872.85, 881.55, 890.7, 900.3, 911.4, 921.15]
                yData = [966.3, 970.2, 975.15, 980.4, 985.8, 992.1, 998.7, 1006.2, 1013.7, 1020.75, 1028.85, 1035.9, 1044.15, 1050.15, 1057.2, 1064.25, 1072.05, 1078.2, 1084.95, 1090.8, 1098.3, 1104.75, 1111.05, 1116.75, 1123.65, 1130.4, 1135.8, 1141.95, 1148.25, 1153.5, 1159.5, 1164.45, 1168.8, 1173.75, 1179, 1182.9, 1187.85, 1191.3, 1195.05, 1199.1, 1202.4, 1205.25, 1207.8, 1208.85, 1210.65, 1212.3, 1213.2, 1214.1, 1214.4, 1215.3, 1216.8, 1217.1, 1217.7, 1217.85, 1218, 1217.1, 1216.35, 1215, 1214.4, 1213.2, 1212.3, 1211.1, 1209.6, 1207.35, 1204.2, 1200.6, 1196.25, 1191.9, 1187.25, 1182.9, 1178.1, 1174.05, 1168.8, 1162.95, 1158.3, 1154.1, 1148.55, 1143.6, 1139.25, 1135.5, 1130.85, 1125.9, 1118.7, 1113.6, 1109.1, 1103.7, 1098.15, 1094.7, 1089.6, 1086.45, 1082.4, 1077.3, 1073.1, 1068.9, 1064.1, 1059.45, 1055.1, 1050.45, 1046.4, 1042.5, 1038.3, 1033.95, 1029.3, 1024.35, 1021.35, 1017, 1012.95, 1008.9, 1005.45, 1002.75, 1000.65, 997.35, 995.1, 993.15, 991.5, 988.95, 986.85, 983.4, 980.7, 978.6, 976.2, 973.65, 972.45, 970.5, 969.45, 968.25, 966.3, 963.75, 962.1, 959.7, 958.5, 957.3, 955.95, 955.35, 954.15, 955.05, 954.3, 954.45, 953.85, 954.6, 954.15, 955.95, 955.8, 957.15, 957.75, 959.85, 961.65, 964.5, 966.6, 969.3, 972.9, 977.25, 981.15, 984.6, 989.4, 993.9, 998.1, 1003.2, 1008.3, 1013.4, 1019.55, 1025.4, 1033.5, 1041.15, 1049.25, 1057.35, 1066.95, 1075.2, 1085.4, 1093.05, 1102.2, 1110.15]
                zData = [-740.85, -740.55, -741.15, -741.75, -742.5, -743.4, -745.05, -748.5, -751.95, -756.15, -759.45, -764.1, -769.8, -775.8, -781.2, -787.35, -794.4, -801.3, -809.85, -817.8, -826.05, -834.45, -843.3, -852.15, -861.45, -870, -878.85, -886.8, -896.1, -904.95, -912.45, -920.85, -928.2, -935.1, -943.05, -949.8, -955.95, -964.2, -970.8, -978.45, -985.8, -994.05, -1000.65, -1007.7, -1014.45, -1021.65, -1028.1, -1035.15, -1040.25, -1045.8, -1051.65, -1056.6, -1060.65, -1064.55, -1067.55, -1071, -1074.75, -1078.05, -1081.2, -1084.35, -1087.2, -1090.5, -1093.05, -1095.3, -1097.7, -1098.45, -1099.5, -1099.8, -1100.25, -1100.4, -1099.5, -1097.55, -1095.9, -1094.1, -1093.05, -1090.65, -1088.25, -1085.25, -1082.25, -1078.95, -1076.1, -1071.45, -1066.5, -1062.6, -1058.55, -1053.6, -1049.1, -1044.6, -1041.15, -1036.5, -1031.85, -1027.8, -1024.65, -1019.55, -1014.3, -1008.45, -1004.7, -999.45, -993.9, -988.5, -982.95, -977.1, -970.05, -963.15, -956.55, -950.4, -943.5, -938.1, -931.65, -926.7, -921.75, -915.45, -910.05, -905.4, -900.9, -897.45, -893.7, -889.2, -885.6, -881.55, -877.5, -873.3, -868.8, -864.45, -859.65, -854.7, -849.15, -843.9, -837.45, -831.6, -825.6, -819.3, -813.3, -807.9, -801.6, -797.1, -792.6, -787.8, -784.05, -780, -776.1, -773.25, -770.25, -767.25, -764.1, -761.25, -759, -757.5, -756.3, -754.8, -753.15, -753.45, -753.45, -752.85, -753.45, -754.35, -755.1, -756.3, -757.35, -759.45, -763.65, -766.8, -770.1, -775.8, -781.95, -789.15, -797.4, -805.35, -815.1, -824.7, -833.85, -843.45]
                uData = [1536.45, 1519.35, 1364.7, 1144.5, 1008.75, 1023.6, 1181.85, 1400.7, 1531.5]
                vData = [1518.6, 1443.75, 1313.25, 1192.95, 1150.65, 1180.65, 1316.4, 1439.85, 1469.1]
                break

            case "angled": // mounted at strange angle with spin axis components []: spin-axis 45 degrees to X, Y & Z
                scanTimes = [23173, 23193, 23213, 23233, 23253, 23273, 23293, 23313, 23333, 23353, 23373, 23393, 23413, 23433, 23453, 23473, 23493, 23513, 23535, 23557, 23577, 23597, 23617, 23637, 23657, 23677, 23697, 23717, 23737, 23757, 23777, 23797, 23817, 23837, 23857, 23877, 23897, 23917, 23937, 23957, 23977, 23997, 24017, 24037, 24057, 24077, 24097, 24117, 24137, 24157, 24177, 24197, 24217, 24237, 24257, 24277, 24297, 24317, 24337, 24357, 24377, 24397, 24417, 24437, 24457, 24481, 24501, 24521, 24541, 24561, 24581, 24601, 24621, 24641, 24661, 24681, 24701, 24721, 24741, 24761, 24781, 24801, 24821, 24841, 24861, 24881, 24901, 24921, 24941, 24961, 24981, 25001, 25021, 25041, 25061, 25081, 25101, 25121, 25141, 25165, 25185, 25205, 25225, 25245, 25265, 25285, 25305, 25329, 25349, 25369, 25389, 25409, 25429, 25449, 25469, 25489, 25509, 25532, 25553, 25573, 25593, 25613, 25633, 25653, 25673, 25693, 25713, 25733, 25753, 25773, 25793, 25813, 25833, 25853, 25873, 25893, 25913, 25933, 25953, 25976, 25997, 26017, 26037, 26057, 26077, 26097, 26117, 26137, 26157, 26177, 26197, 26217, 26237, 26257, 26277, 26300, 26321, 26341, 26361, 26381, 26401, 26421, 26441, 26461, 26481, 26501, 26521, 26544, 26565, 26585, 26605, 26625, 26645, 26665, 26685, 26705, 26729, 26753, 26773, 26793, 26813, 26833, 26853, 26873, 26893, 26913, 26933, 26953, 26973, 26993, 27013, 27033, 27053, 27073, 27093, 27113, 27133, 27153, 27181, 27201, 27221, 27241, 27261, 27281, 27301, 27321, 27341, 27361, 27381, 27401, 27421, 27441, 27461, 27481, 27501, 27521, 27541, 27561, 27581, 27601, 27621, 27641, 27661, 27681, 27701, 27721, 27741, 27761, 27781, 27801, 27821, 27841, 27861, 27881, 27901, 27921, 27941, 27961, 27981, 28001, 28021, 28041]
                xData = [-202.95, -205.35, -207.45, -210.3, -211.35, -212.55, -214.2, -216.9, -219.15, -222.45, -224.55, -226.65, -229.2, -231, -231.6, -232.5, -232.8, -233.4, -234.15, -234.15, -235.05, -235.05, -234.6, -234.3, -233.25, -230.85, -229.2, -226.65, -224.55, -222.15, -219.15, -216.15, -214.05, -211.2, -207.45, -203.25, -198.3, -193.2, -189, -184.2, -178.35, -174.3, -169.95, -165.75, -161.55, -157.5, -153, -148.95, -143.25, -137.25, -133.2, -128.7, -122.25, -117.3, -111.45, -106.2, -102.15, -96.75, -90.6, -86.55, -81.15, -77.25, -73.05, -67.95, -63.6, -61.5, -57.15, -54.3, -50.4, -46.65, -44.7, -40.8, -36.15, -34.05, -31.2, -29.4, -27.45, -24.6, -23.55, -22.95, -21.9, -20.1, -18.6, -18.15, -16.95, -16.35, -16.2, -16.05, -16.65, -16.95, -17.7, -20.25, -20.85, -22.05, -24, -24.6, -27.15, -28.35, -29.85, -32.7, -35.25, -36.9, -41.4, -43.95, -47.7, -51, -54.75, -58.8, -64.65, -69.75, -76.05, -82.65, -89.1, -95.25, -101.4, -106.65, -111.15, -115.05, -120.75, -125.55, -130.95, -134.7, -139.35, -144.3, -148.8, -151.8, -154.95, -159.6, -164.25, -167.7, -171.3, -176.1, -179.85, -184.05, -186.9, -190.05, -194.1, -197.4, -199.5, -202.5, -206.25, -208.95, -211.8, -213.3, -216.15, -219, -220.95, -222.6, -224.7, -226.05, -228.6, -229.8, -231.3, -232.5, -232.65, -232.5, -232.5, -232.5, -231.6, -231.45, -231.6, -231.3, -229.65, -228.6, -225.45, -223.8, -220.2, -216.75, -212.85, -210.6, -206.85, -204, -199.8, -196.5, -191.85, -188.4, -183.45, -179.1, -173.85, -169.5, -164.4, -160.65, -155.55, -151.05, -146.1, -141.75, -138, -133.95, -129.6, -126, -121.65, -118.05, -113.85, -109.05, -105, -100.05, -95.55, -90.45, -85.65, -81.9, -77.7, -73.05, -68.4, -63.9, -60.3, -57.45, -52.95, -49.2, -45.3, -43.2, -40.2, -37.95, -34.65, -31.8, -30.15, -27.45, -24.75, -22.05, -20.25, -17.7, -17.4, -15.15, -15.6, -14.7, -13.8, -13.5, -14.4, -14.25, -15, -15.45, -16.65, -18.3, -19.8, -21.15, -21.75, -22.5, -23.85, -25.05, -27.45, -29.25, -30.75, -34.35]
                yData = [230.4, 233.55, 236.55, 239.55, 242.1, 245.25, 247.8, 251.4, 254.4, 258, 262.2, 267.45, 272.25, 277.65, 282.3, 288, 292.65, 298.35, 301.8, 307.05, 311.85, 316.8, 321.3, 325.95, 329.85, 334.2, 338.25, 342.3, 347.25, 350.85, 355.2, 358.95, 364.5, 367.2, 370.8, 372.45, 376.8, 379.95, 382.35, 383.85, 387.15, 388.8, 391.65, 391.8, 392.4, 392.55, 393.45, 393.45, 394.8, 394.35, 394.05, 393.9, 393.15, 392.4, 391.35, 388.2, 386.7, 385.35, 383.1, 381.75, 378.75, 376.5, 374.4, 371.25, 369.3, 366, 363, 359.85, 356.4, 353.85, 350.25, 344.55, 341.4, 336.45, 331.8, 326.85, 321.75, 317.55, 313.8, 309, 304.5, 300.3, 295.8, 291.15, 286.05, 281.4, 275.7, 271.05, 265.8, 261.75, 256.95, 252.3, 248.4, 244.5, 241.2, 237.45, 232.95, 228.75, 225.45, 221.55, 218.55, 215.25, 212.1, 208.8, 206.4, 203.1, 200.1, 197.7, 194.55, 193.05, 191.55, 189.9, 188.85, 187.5, 187.2, 187.65, 186.75, 187.65, 188.55, 189.3, 190.2, 190.5, 191.4, 192.9, 194.4, 195.15, 196.65, 198.6, 200.4, 202.95, 205.05, 206.55, 209.7, 212.4, 215.7, 219.3, 220.2, 223.8, 226.5, 228.3, 231.6, 234, 235.8, 240.6, 243.6, 247.8, 253.05, 257.1, 262.5, 267.9, 272.25, 278.4, 283.5, 287.7, 293.7, 298.95, 304.8, 310.95, 316.8, 322.95, 329.25, 333.3, 337.35, 342.3, 346.05, 348.9, 352.65, 357.3, 361.5, 365.85, 369.3, 372.45, 375.75, 378.45, 380.25, 382.95, 385.65, 387.3, 390, 391.65, 393.45, 394.65, 396.3, 397.2, 398.25, 398.1, 398.4, 398.1, 397.8, 397.35, 395.55, 394.05, 393.6, 393.15, 390.9, 388.65, 386.1, 384.45, 382.65, 379.95, 377.25, 375.15, 373.2, 369.6, 366.15, 363.3, 360.15, 356.25, 351.75, 348.15, 345.45, 342.15, 337.2, 333.3, 330, 326.85, 322.05, 317.25, 312.3, 308.55, 303.3, 297.9, 292.95, 289.2, 284.55, 280.65, 275.55, 272.25, 267.3, 264.3, 259.05, 255.15, 250.8, 247.35, 242.85, 239.7, 235.35, 232.5, 229.2, 226.05, 223.2, 220.2]
                zData = [456.6, 457.2, 456.45, 456.15, 456.75, 456.6, 455.55, 455.85, 454.35, 453.6, 450.6, 448.5, 445.35, 443.85, 441.9, 439.05, 436.5, 436.2, 433.65, 431.7, 428.85, 426, 423.9, 421.2, 417, 413.25, 410.1, 407.7, 405.45, 402.6, 399.3, 396.9, 394.05, 391.65, 387.9, 383.85, 382.2, 379.35, 376.5, 373.8, 371.25, 369, 366.9, 363.15, 359.7, 357.3, 355.35, 352.95, 350.85, 348.3, 347.25, 346.65, 344.25, 342.45, 340.35, 338.55, 338.25, 336.45, 334.5, 333.45, 331.2, 330.45, 330.9, 329.85, 329.55, 330.15, 330.75, 332.25, 332.25, 331.65, 332.1, 332.25, 332.85, 333.6, 335.1, 336.15, 337.5, 338.7, 340.05, 342.45, 344.7, 347.85, 350.85, 353.4, 356.55, 359.4, 360.3, 362.4, 363.45, 366.15, 368.7, 371.7, 374.25, 377.25, 379.8, 382.35, 385.05, 388.05, 389.85, 393.3, 397.95, 401.7, 404.25, 406.35, 409.2, 412.65, 415.5, 416.7, 419.25, 423, 427.05, 430.5, 433.5, 436.2, 439.65, 442.5, 445.05, 446.1, 446.7, 447.9, 449.7, 451.2, 451.8, 453.45, 455.1, 455.4, 455.85, 455.4, 456.3, 456.3, 455.1, 456.6, 458.85, 459.45, 459.45, 459.15, 459.9, 460.2, 458.4, 457.95, 457.05, 457.2, 455.85, 454.5, 454.05, 452.55, 450.6, 450.45, 448.8, 448.05, 447, 444.3, 443.55, 442.35, 439.8, 437.25, 434.1, 431.1, 428.7, 425.55, 422.25, 418.8, 415.5, 411.9, 408.3, 404.7, 401.55, 398.7, 396, 394.35, 392.1, 389.25, 387.6, 384.6, 381.3, 377.4, 373.35, 370.8, 367.8, 363.3, 360.9, 357.9, 355.95, 354.45, 351.9, 350.1, 349.05, 346.65, 344.55, 343.2, 341.7, 340.2, 338.4, 337.95, 336.15, 335.7, 335.1, 334.05, 332.85, 332.7, 331.35, 332.25, 331.2, 329.85, 329.55, 330.6, 330.9, 331.65, 331.95, 333.6, 335.1, 336.3, 336.9, 338.4, 339, 340.2, 340.5, 341.85, 343.2, 344.85, 346.2, 348.45, 349.8, 351.3, 354.15, 355.95, 357.45, 359.7, 362.55, 365.4, 370.05, 371.55, 373.65, 376.2, 378.45, 380.25, 382.35, 383.55, 386.85, 389.55, 391.8, 393.9]
                uData = [-163.9, -219.44, -220.54, -166.76, -85.89, -24.94, -33.13, -85.38, -165.86]
                vData = [192.23, 262.34, 351, 409.81, 394.16, 319.22, 233.76, 182.8, 193.04]
                
                break

            default: // tetrahedral: spin-axis 45 degrees to X, Y & Z
                scanTimes = [23173, 23193, 23213, 23233, 23253, 23273, 23293, 23313, 23333, 23353, 23373, 23393, 23413, 23433, 23453, 23473, 23493, 23513, 23535, 23557, 23577, 23597, 23617, 23637, 23657, 23677, 23697, 23717, 23737, 23757, 23777, 23797, 23817, 23837, 23857, 23877, 23897, 23917, 23937, 23957, 23977, 23997, 24017, 24037, 24057, 24077, 24097, 24117, 24137, 24157, 24177, 24197, 24217, 24237, 24257, 24277, 24297, 24317, 24337, 24357, 24377, 24397, 24417, 24437, 24457, 24481, 24501, 24521, 24541, 24561, 24581, 24601, 24621, 24641, 24661, 24681, 24701, 24721, 24741, 24761, 24781, 24801, 24821, 24841, 24861, 24881, 24901, 24921, 24941, 24961, 24981, 25001, 25021, 25041, 25061, 25081, 25101, 25121, 25141, 25165, 25185, 25205, 25225, 25245, 25265, 25285, 25305, 25329, 25349, 25369, 25389, 25409, 25429, 25449, 25469, 25489, 25509, 25532, 25553, 25573, 25593, 25613, 25633, 25653, 25673, 25693, 25713, 25733, 25753, 25773, 25793, 25813, 25833, 25853, 25873, 25893, 25913, 25933, 25953, 25976, 25997, 26017, 26037, 26057, 26077, 26097, 26117, 26137, 26157, 26177, 26197, 26217, 26237, 26257, 26277, 26300, 26321, 26341, 26361, 26381, 26401, 26421, 26441, 26461, 26481, 26501, 26521, 26544, 26565, 26585, 26605, 26625, 26645, 26665, 26685, 26705, 26729, 26753, 26773, 26793, 26813, 26833, 26853, 26873, 26893, 26913, 26933, 26953, 26973, 26993, 27013, 27033, 27053, 27073, 27093, 27113, 27133, 27153, 27181, 27201, 27221, 27241, 27261, 27281, 27301, 27321, 27341, 27361, 27381, 27401, 27421, 27441, 27461, 27481, 27501, 27521, 27541, 27561, 27581, 27601, 27621, 27641, 27661, 27681, 27701, 27721, 27741, 27761, 27781, 27801, 27821, 27841, 27861, 27881, 27901, 27921, 27941, 27961, 27981, 28001, 28021, 28041]
                xData = [-160.80, -184.05, -191.85, -198.90, -207.30, -213.60, -221.70, -227.85, -232.35, -237.75, -243.00, -248.25, -252.45, -255.90, -260.85, -265.05, -268.65, -272.70, -275.10, -279.00, -281.70, -283.20, -285.75, -287.55, -288.45, -289.95, -290.70, -291.15, -291.75, -291.15, -291.00, -289.95, -289.35, -288.00, -287.55, -286.05, -285.00, -283.05, -282.00, -279.90, -277.65, -274.50, -270.30, -266.40, -261.75, -256.05, -251.25, -247.05, -241.95, -237.30, -232.35, -227.40, -223.95, -218.55, -212.85, -207.45, -203.70, -197.85, -192.30, -186.90, -181.65, -175.65, -170.25, -163.05, -156.00, -150.60, -144.30, -138.15, -132.30, -126.45, -121.50, -116.55, -111.15, -105.15, -100.20, -95.85, -90.75, -84.60, -80.10, -75.30, -72.00, -68.10, -63.60, -59.55, -57.15, -54.90, -52.95, -50.55, -48.30, -47.40, -47.25, -47.85, -48.15, -48.90, -50.10, -53.10, -55.80, -58.95, -61.20, -64.35, -67.95, -72.15, -76.05, -80.40, -83.40, -87.75, -91.65, -94.20, -97.65, -102.00, -104.85, -109.80, -112.65, -117.90, -124.95, -131.10, -135.90, -142.35, -147.30, -153.90, -159.75, -164.70, -169.95, -174.60, -180.00, -186.75, -192.90, -197.85, -202.80, -206.40, -211.35, -214.80, -219.30, -223.20, -226.35, -229.80, -234.45, -238.20, -242.10, -244.80, -248.40, -253.20, -256.80, -260.55, -264.30, -267.30, -271.80, -274.65, -277.20, -280.50, -282.45, -285.45, -288.60, -289.50, -290.55, -291.45, -292.80, -292.35, -292.05, -292.05, -291.45, -291.45, -290.55, -287.70, -287.70, -286.80, -284.85, -281.70, -279.30, -275.85, -273.75, -269.55]
                yData = [2.55, 2.55, 3.15, 4.5, 4.05, 4.35, 5.25, 4.8, 3.75, 4.2, 3.9, 4.2, 3.9, 3.6, 3.3, 3.9, 2.4, 2.25, 1.5, 0.6, -1.05, -2.85, -4.35, -5.4, -8.1, -8.7, -10.05, -10.65, -11.4, -13.05, -14.55, -15, -18.3, -20.7, -22.95, -25.2, -27.15, -30.45, -34.05, -36.6, -38.25, -40.8, -43.2, -47.1, -48.6, -51, -53.55, -57.15, -61.05, -64.8, -67.35, -71.85, -76.65, -81.45, -85.95, -90.45, -94.35, -98.25, -102.3, -105.9, -109.35, -112.05, -113.55, -115.2, -115.8, -115.2, -113.25, -112.2, -111.45, -111.75, -113.1, -114.45, -116.4, -120.3, -123.45, -126.75, -129.6, -131.55, -133.95, -135.6, -136.5, -136.95, -138.45, -139.5, -140.25, -141.9, -143.25, -144.15, -144.6, -144.6, -144.75, -145.65, -145.8, -146.1, -146.4, -147.9, -148.35, -148.8, -148.95, -148.8, -149.85, -150.75, -151.5, -152.4, -152.55, -152.4, -153.3, -153, -153.15, -151.65, -150.9, -152.25, -152.4, -152.25, -152.1, -151.35, -151.35, -151.2, -149.1, -148.2, -146.55, -145.8, -145.35, -146.1, -144.9, -144.9, -144.75, -144.45, -144, -142.95, -140.85, -141, -139.65, -139.35, -139.05, -137.7, -136.65, -136.35, -134.25, -133.95, -131.55, -130.05, -129.9, -129.3, -127.5, -127.05, -125.55, -125.25, -124.5, -123, -121.95, -120.75, -119.4, -119.25, -118.35, -117.75, -115.65, -113.85, -111.6, -110.7, -107.85, -106.35, -103.2, -101.55, -99.3, -98.25, -95.7, -94.35, -91.65, -91.05, -89.4, -89.1, -87.9, -87.45, -85.65, -85.2, -83.55, -82.2, -80.25, -78.15, -76.2, -74.4, -70.65, -68.1, -66, -61.95, -59.4, -54.9, -51.75, -50.1, -47.1, -44.1, -42.6, -41.1, -39.45, -37.2, -34.8, -32.7, -30.9, -28.65, -25.2, -23.4, -22.2, -21.45, -20.1, -18, -17.1, -17.4, -16.5, -15.3, -13.8, -13.05, -12.6, -10.8, -8.85, -7.65, -5.1, -2.85, -0.6, 0.75, 2.25, 3.6, 4.95, 5.7, 6.15, 7.05, 7.65, 8.25, 9, 9.6, 10.05, 11.4, 11.7, 13.05, 13.35, 15, 15.6, 15.75, 14.7, 14.55, 14.85, 14.7, 13.65, 13.65, 13.8, 14.4, 14.55, 14.1, 14.7, 14.55, 14.55, 15, 15, 14.85, 14.85, 14.1, 13.8, 13.5, 12.3, 11.55, 12.3, 12, 12.45, 12, 11.1, 10.8, 11.1, 10.65, 10.35, 10.65, 10.65, 10.95, 9.9, 9.3, 7.65, 6.6, 5.25, 4.2, 2.4, 1.95, 1.35, 0.9, 0.45, -1.5, -2.85, -3.45, -3.15, -4.95, -6.75, -8.7, -10.2, -11.4, -12.75, -15.6, -16.05, -16.65, -16.8, -17.55, -18.45, -19.65, -19.8, -21.45, -22.65, -25.05, -26.85, -28.35, -30.45, -32.4, -33.15, -35.55, -37.05, -38.25, -39.6, -39.6, -41.7, -45.3, -46.8, -48.9, -51.15, -54.15, -58.35, -61.05, -63.3, -66.3, -69, -71.7, -75, -77.85, -81, -84.6, -87.6, -91.5, -95.55, -97.8, -100.95, -103.8, -106.5, -108.75, -110.85, -112.65, -115.35, -117, -119.4, -121.65, -124.2, -126, -127.5, -129.45, -131.55, -132.9, -133.8, -135, -136.95, -139.2, -141.3, -142.8, -144.45, -147.15, -149.85, -151.05, -152.25, -152.85, -153.3, -154.8, -154.05, -154.05, -154.8, -154.65, -153.75, -153.15, -152.4, -152.55, -151.95, -150.6, -150.3, -149.55, -149.7, -149.4, -148.2, -147.3, -147.3, -146.55, -146.25, -144.75, -143.1, -142.05, -140.55, -138, -135.9, -135, -133.95, -133.05, -131.1, -130.8, -130.35, -129.9, -129.3, -128.4, -127.35, -127.35, -125.7, -124.65, -123.75, -122.55, -122.7, -121.65, -120.75, -119.7, -119.1, -118.35, -116.85, -114.15, -112.95, -111.3, -110.4, -108.75, -106.8, -105, -103.65, -101.4, -98.4, -94.8, -91.05, -87.9, -84.6, -81.75, -78.3, -76.35, -74.7, -73.8, -71.7, -69.75, -68.55, -67.35, -65.85, -63.6, -62.7, -61.95, -61.2, -59.1, -57.9, -56.4, -54.45, -51.3, -48.6, -45.3, -42.6, -39.3, -35.85, -34.35, -32.55, -31.65, -30.75, -28.65, -28.2, -27.3, -25.65, -23.55, -20.85, -19.8, -18.45, -15.75, -14.4, -12, -9, -6.15, -2.4, -0.15, 1.65, 3.75, 5.55, 5.55, 6.45, 7.05, 8.4, 9.75, 10.95, 11.7, 13.8, 14.25, 14.55, 13.95, 14.25, 14.1, 13.95, 13.5, 14.4, 14.1, 14.4, 13.35, 13.35, 13.05, 12.15, 10.95, 11.25, 11.7, 11.7, 12.15, 11.85, 12.75, 12.9, 13.05, 12.9, 12.45, 11.4, 10.8, 10.5, 8.85, 7.5, 6.15, 4.95, 4.35, 3.15, 1.05, 0.6, -1.05, -2.7, -3, -5.25, -6.3]
                zData = [379.2, 377.7, 374.7, 369.15, 366.15, 362.85, 359.25, 356.85, 353.4, 349.8, 346.8, 343.65, 338.85, 334.8, 330.6, 325.5, 321, 316.8, 312.9, 310.2, 307.95, 303, 300.9, 297.45, 294.9, 291.6, 289.65, 286.65, 285.15, 282.3, 280.8, 277.95, 274.8, 269.85, 265.5, 262.5, 258.6, 254.4, 252.45, 249.6, 247.95, 247.05, 245.85, 244.65, 243.9, 242.55, 242.1, 241.35, 240.45, 239.55, 237.45, 236.55, 234.6, 234, 233.7, 233.85, 233.7, 235.65, 236.1, 238.35, 238.35, 240.15, 240.75, 240.6, 240.75, 241.65, 240.9, 241.2, 241.2, 241.8, 244.2, 245.55, 247.65, 249.9, 253.65, 255.15, 258.45, 260.25, 262.65, 263.7, 265.8, 266.25, 267.9, 267.9, 269.1, 270.45, 272.55, 273, 274.95, 275.7, 277.8, 279.3, 281.4, 283.65, 287.4, 290.1, 294.9, 298.2, 302.25, 305.85, 309, 312.6, 316.2, 319.35, 322.65, 326.4, 328.8, 332.25, 333.9, 337.05, 338.25, 339.75, 340.35, 342, 343.5, 345.75, 346.95, 350.1, 352.8, 356.1, 358.05, 359.55, 362.25, 364.65, 366.3, 368.7, 369.45, 371.85, 373.8, 374.85, 375.75, 377.85, 378.45, 380.25, 382.8, 385.35, 387.45, 388.65, 390.3, 391.2, 392.4, 392.55, 393.75, 394.65, 397.5, 397.95, 400.95, 402.15, 403.8, 403.5, 405.3, 406.35, 408.3, 409.2, 410.85, 412.05, 415.35, 416.85, 418.05, 419.1, 421.5, 423.75, 426.3, 427.05, 428.1, 429.75, 431.4, 431.1, 432.3, 431.4, 433.05, 433.8, 435.15, 435, 435.3, 435.45, 437.25, 437.1, 437.1, 436.35, 436.8, 437.7, 437.7, 437.7, 437.85, 438.3, 438, 438.15, 437.4, 436.65, 436.05, 435.3, 433.65, 433.65, 432.6, 432.15, 431.25, 430.35, 429.6, 429.6, 427.05, 425.85, 424.35, 422.7, 421.5, 420.6, 418.65, 418.8, 417.45, 416.55, 416.55, 414.45, 412.8, 411.3, 408.3, 406.35, 403.8, 400.95, 398.55, 396, 393.6, 390.6, 388.8, 386.7, 383.7, 381.6, 378.3, 376.2, 374.85, 371.4, 368.55, 365.7, 363, 360.6, 357.75, 355.5, 353.25, 351.75, 350.55, 349.05, 348.75, 347.25, 345.9, 345.45, 345, 345.15, 344.1, 342.15, 342.45, 342.15, 342, 340.35, 338.55, 337.8, 337.2, 335.25, 333.75, 331.65, 329.25, 327, 324.75, 321.75, 319.35, 317.55, 314.85, 312.6, 310.5, 308.4, 306.6, 305.25, 302.1, 301.5, 300.9, 299.55, 298.35, 297.3, 296.1, 294.9, 291.45, 288.9, 286.5, 285.15, 283.8, 280.65, 278.7, 277.5, 276.6, 274.35, 272.7, 269.85, 268.8, 267.45, 265.65, 263.55, 262.8, 260.7, 260.85, 259.65, 258.15, 257.4, 256.35, 254.25, 252.45, 249.75, 249.3, 248.4, 247.5, 246.45, 245.85, 244.8, 244.2, 241.95, 240.6, 239.7, 238.2, 237.9, 237, 235.5, 234.45, 233.25, 232.2, 232.05, 231.45, 231.6, 232.2, 231.45, 231.75, 231.45, 229.95, 230.4, 228.9, 228.9, 231, 232.35, 232.2, 234.75, 235.65, 238.5, 238.95, 240.3, 241.35, 244.05, 246.45, 248.4, 249.6, 252.3, 252.6, 253.65, 254.85, 254.7, 256.35, 259.05, 261.6, 266.1, 269.85, 274.95, 279.9, 284.7, 289.5, 295.35, 299.7, 305.4, 310.2, 316.05, 322.5, 327.75, 333, 337.65, 340.8, 344.7, 347.4, 349.05, 351.45, 352.95, 354.3, 356.7, 358.2, 360.75, 361.95, 363.9, 365.55, 368.25, 370.05, 372.75, 375.45, 379.05, 382.5, 385.05, 387.3, 390.45, 393.15, 394.35, 396.45, 396.9, 398.4, 400.35, 401.85, 401.85, 403.65, 403.65, 404.1, 404.7, 404.55, 405.45, 406.5, 406.95, 407.55, 409.35, 410.55, 412.8, 413.25, 414.15, 416.1, 416.7, 417.45, 418.95, 420, 421.2, 423.3, 424.35, 427.8, 430.35, 431.7, 434.25, 436.8, 438.45, 439.2, 440.1, 440.4, 441.45, 440.25, 439.95, 439.5, 440.1, 440.25, 440.7, 440.85, 442.35, 442.35, 442.35, 442.5, 441.6, 441.3, 439.8, 438.6, 437.4, 436.2, 434.1, 433.05, 430.5, 429.9, 428.25, 427.05, 426.9, 425.7, 423.9, 423.45, 421.05, 418.5, 416.1, 413.7, 410.85, 408.45, 405.6, 403.05, 400.65, 398.25, 394.65, 391.8, 388.5, 384.6, 381.6, 378.6, 375.45, 371.1, 368.4, 366, 363.6, 361.35, 359.25, 356.7, 354.6, 353.1, 350.85, 350.1, 348.75, 347.4, 346.5, 347.1, 345.45, 344.55, 343.05, 340.65, 339, 337.2, 335.7, 334.05, 332.25, 328.2, 325.95, 323.55, 321.6, 318, 314.55, 311.85, 310.2, 307.95, 305.55, 302.25, 300.3, 298.5, 295.8, 294.3, 291.9, 288.6, 286.05, 283.2, 281.1, 277.35, 274.35]
                uData = [-829.05, -829.80, -828.00, -837.45, -836.40, -856.65, -862.20, -860.85, -860.10, -868.95, -862.50, -861.15]
                vData = [-367.35, -360.60, -181.65, -435.30, -251.85, -283.35, -288.00, -284.40, -317.25, -341.40, -237.30, -209.55]
                break
        }

        // transpose the three arrays into array of triples
        for (let n = 0; n < scanTimes.length; n++) {
            let xyz = []
            xyz.push(xData[n])
            xyz.push(yData[n])
            xyz.push(zData[n])
            scanData.push(xyz)
        }
    }

    // helpful for logging:
    function round2(v: number): number {
        return (Math.round(100 * v) / 100)
    }

}
