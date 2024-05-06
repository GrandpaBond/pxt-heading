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
        cloneMe(): Arrow {
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
                    if (Math.abs(deviate) < (Math.PI / 2)) {
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
                    period = (endTime - startTime) / turns
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
                    majors.splice(i, 1)  // disqualified!
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
    let capturing = false // data capturing mode
    let debugging = false  // test mode flag
    let dataset: string = "NONE" // test dataset to use
    let testData: number[][] = [] //[X,Y,Z] magnetometer values for test cases
    let test = 0 // global selector for test-cases
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
        if ((views[View.XY].period == -1)
            && (views[View.YZ].period == -1)
            && (views[View.ZX].period == -1)) {
            return -3 // "NOT ENOUGH SCAN ROTATION"
        }

        // Choose the "roundest" Ellipse  --the one with lowest eccentricity.
        bestView = View.XY
        if (views[View.YZ].eccentricity < views[bestView].eccentricity) bestView = View.YZ
        if (views[View.ZX].eccentricity < views[bestView].eccentricity) bestView = View.ZX

        // periodicity is unreliable in the best View: average just the other two Views' measurements
        period = (views[0].period + views[1].period + views[2].period - views[bestView].period) / 2
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

    /* Although eventually we'd only need [uDim, vDim], we sum and log all three dims.
       This will allow us, while testing, to override automatic choice of bestView
    */

    function takeSingleReading(): number {
        let uRaw = 0
        let vRaw = 0
        let u = 0
        let v = 0
        let uNew = 0
        let vNew = 0
        let uFix = 0
        let vFix = 0

        let reading = 0
        if (debugging) { // just choose the next test-data value
            uRaw = testData[test][uDim]
            vRaw = testData[test][vDim]
            test = (test + 1) % testData.length
        } else {
            let xyz = [0, 0, 0]

            // get a new sample as the sum of seven readings
            xyz[0] = input.magneticForce(0)
            xyz[1] = input.magneticForce(1)
            xyz[2] = input.magneticForce(2)

            for (let i = 0; i < 6; i++) {
                basic.pause(5)
                xyz[0] += input.magneticForce(0)
                xyz[1] += input.magneticForce(1)
                xyz[2] += input.magneticForce(2)
            }

            if (capturing) {
                datalogger.log(
                    datalogger.createCV("index", test),
                    datalogger.createCV("x", round2(xyz[0])),
                    datalogger.createCV("y", round2(xyz[1])),
                    datalogger.createCV("z", round2(xyz[2])))
            }

            testData.push(xyz)
            // use the test sample we just captured (in our current bestView)
            uRaw = testData[test][uDim]
            vRaw = testData[test][vDim]

            test++ // clock another test sample
        }
        // normalise the point w.r.t our Ellipse origin
        u = uRaw - uOff
        v = vRaw - vOff

        // Unless this Ellipse.isCircular, any {u,v} reading will be foreshortened in this View, and
        // must stretched along the Ellipse minor-axis, to place it correctly on the Spin-Circle.
        if (needsFixing) {
            // First rotate CLOCKWISE by theta (so aligning the Ellipse minor-axis angle with the V-axis)
            uNew = u * cosTheta + v * sinTheta
            vNew = v * cosTheta - u * sinTheta
            uFix = uNew
            // Now scale up vertically, re-balancing the axes to make the Ellipse circular
            vFix = vNew * scale
            // get the adjusted angle for this corrected {u,v}
            reading = Math.atan2(vFix, uFix)
            // finally, undo the rotation by theta
            reading += theta
        } else {
            reading = Math.atan2(v, u)
        }

        if (debugging) {
            datalogger.log(
                datalogger.createCV("u", round2(u)),
                datalogger.createCV("v", round2(v)),
                datalogger.createCV("uNew", round2(uNew)),
                datalogger.createCV("vNew", round2(vNew)),
                datalogger.createCV("uFix", round2(uFix)),
                datalogger.createCV("vFix", round2(vFix)),
                datalogger.createCV("reading", Math.round(reading)),
                datalogger.createCV("bearing", Math.round(asBearing(reading)))
            )
        }
        return reading
    }

    // helpful for logging...
    function round2(v: number): number {
        return (Math.round(100 * v) / 100)
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
    // that were captured from live runs using the datalogger.
    function simulateScan(dataset: string) {
        let xData: number[] = []
        let yData: number[] = []
        let zData: number[] = []
        let xTest: number[] = []
        let yTest: number[] = []
        let zTest: number[] = []
        switch (dataset) {

            case "angled": // mounted at strange angle with spin axis components []: spin-axis 45 degrees to X, Y & Z
                scanTimes = []
                xData = []
                yData = []
                zData = []
                xTest = []
                yTest = []
                zTest = []
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
        // do the same fo the test cases
        for (let n = 0; n < xTest.length; n++) {
            let xyz = []
            xyz.push(xTest[n])
            xyz.push(yTest[n])
            xyz.push(zTest[n])
            testData.push(xyz)
        }
    }
}
