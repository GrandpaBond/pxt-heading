
 // TODO? No use is yet made of the accelerometer. Although seemingly helpful to compensate for
 // static tilt, on a moving buggy dynamic sideways accelerations confound the measurement of "down",
 // so applying tilt-compensation could actually make compass-heading accuracy worse!


// OPERATIONAL MODES (controlling data-logging)
enum Mode {
    Normal, // Normal usage, mounted in a buggy
    Capture, // Acquire a new test dataset, using a special rotating jig
    Debug, // Test & debug (NOTE: named sets of test-dataset are hard-coded below)
    Trace // Collect full trace (= Capture + Debug) 
}

//% color=#6080e0 weight=40 icon="\uf14e" block="Heading" 
/**
 * An extension providing a compass-bearing for a buggy located anywhere on the globe 
 * (except at the magnetic poles!), with any arbitrary mounting orientation for its microbit.
 * 
 * See the README for a detailed description of the approach, methods and algorithms.
 */
namespace heading {

    // ENUMERATIONS
    enum View { // the axis-pairs that span the three possible magnetometer projection planes
        XY,
        YZ,
        ZX
    }

    // CONSTANTS

    const HalfPi = Math.PI / 2
    const TwoPi = 2 * Math.PI
    const ThreePi = 3 * Math.PI
    const RadianDegrees = 360 / TwoPi
    const EnoughScanTime = 1800 // minimum acceptable scan-time
    const EnoughSamples = 70 // fewest acceptable scan samples
    const TooManySamples = 500 // don't be too greedy with memory!
    const MarginalField = 30 // minimum acceptable field-strength for 7 magnetometer readings
    const Circular = 1.1 // maximum eccentricity to consider an Ellipse as "circular"
    const LongEnough = 0.9 // for major-axis candidates, qualifying fraction of longest radius 

    // SUPPORTING CLASSES

    // An Arrow is an object holding a directed vector {u,v} in both Cartesian and Polar coordinates. 
    // It also carries a time-field, used to timestamp scanned samples.
    // It is used to hold a 2D magnetometer measurement as a re-centred vector.

    class Arrow {
        u: number; // horizontal component
        v: number; // vertical component
        size: number; // polar magnitude of vector
        angle: number; // polar angle (radians anticlockwise from East)
        time: number; // for scan samples, timestamp of when this was collected

        constructor(u: number, v: number, t: number) {
            this.u = u
            this.v = v
            this.size = Math.sqrt((u * u) + (v * v))
            this.angle = 0
            if (this.size > 0) { // defend against zero-divide
                this.angle = Math.atan2(v, u)
            }
            this.time = t
        }
        // when copying...
        cloneMe(): Arrow {
            return new Arrow(this.u, this.v, this.time)
        }
    }

    // An Ellipse is an object holding the characteristics of the (typically) elliptical
    // view formed when projecting the scan Spin-Circle onto a 2-axis View-plane.
    class Ellipse {
        plane: string; // View name (just for debug)
        uDim: number; // horizontal axis of this View
        vDim: number; // vertical axis of this View
        uOff: number; // horizontal offset needed to re-centre this Ellipse along the U-axis
        vOff: number; // vertical offset needed to re-centre this Ellipse along the V-axis

        // calibration characteristics
        majorAxis: Arrow; // direction of the major axis 
        eccentricity: number; // ratio of major-axis to minor-axis magnitudes for this Ellipse
        isCircular: boolean; // flag saying this "Ellipse" View is almost circular, simplifying future handling
        rotationSense: number; // rotation reversal sign = +/-1, reflecting this Ellipse's view of the clockwise scan
        period: number; // scan-rotation time (as twice average time between passing the ellipse major-axis)

        constructor(plane: string, uDim: number, vDim: number, uOff: number, vOff: number) {
            this.plane = plane // (as a DEBUG aid)
            this.uDim = uDim
            this.vDim = vDim
            this.uOff = uOff
            this.vOff = vOff
            this.isCircular = false // until proved otherwise!
            this.period = -1
        }

        // This method analyses a 2-D subset of the scanData for a projected View of the Spin-Circle,
        // creating an Arrow indicating the major-axis direction of the elliptical View.
        // It performs four tasks:
        // 1) By comparing the longest and shortest radii, this method works out the eccentricity of 
        //    the Ellipse, as seen from this View. Also notes whether being viewed "from below".
        // 2) It collects possible candidates for the Ellipse major-axis by looking for local 
        //    radius peaks; candidate values are pushed onto the list of Arrows: this.majors[]
        // 3) It then finds the consensus angle of the axis-candidates (reversing "opposite"
        //    ones, so they all point to the same end of the axis as the first candidate).
        // 4) Works out the average rotation period by clocking each time we pass the first major-axis end. 
        
        analyseView() {
            let majors: Arrow[] = [] // candidate directions for major axis of Ellipse
            let trial = new Arrow(scanData[0][this.uDim], scanData[0][this.vDim], scanTimes[0])
            let longest = trial.size
            let shortest = trial.size
            let spin = 0
            let sizeWas: number
            let angleWas: number
            let step: number = 99999 // marker for "first time round"
            let stepWas: number
            // compare first one with remaining samples
            for (let i = 1; i < scanTimes.length; i++) {
                sizeWas = trial.size
                angleWas = trial.angle
                trial = new Arrow(scanData[i][this.uDim], scanData[i][this.vDim], scanTimes[i])
                longest = Math.max(longest, trial.size)
                shortest = Math.min(shortest, trial.size)
                // accumulate gradual rotation of the projected field-vector...
                let delta = trial.angle - angleWas
                spin += delta
                // fix roll-rounds in either direction
                if (delta > HalfPi) spin -= TwoPi // apparent big positive jump is due to underflow
                if (delta < -HalfPi) spin += TwoPi // apparent big negative jump is due to overflow
                // now collect candidates for the major-axis angle   
                stepWas = step
                step = trial.size - sizeWas // is radius growing or shrinking?
                // (ensure that the first two steps will always match)
                if (stepWas == 99999) stepWas = step
                // look for peaks, where we switch from growing to shrinking
                if ((stepWas > 0) && (step <= 0)) {
                    majors.push(trial.cloneMe()) // copy the major axis we are passing
                }
            }
            // During a clockwise scan, the projection-angle of the field-vector as viewed from above 
            // will appear to rotate anticlockwise (increasing radians) from the perspective of the Buggy. 
            // If we find it rotates clockwise (decreasing radians) then this projection-plane must be 
            // viewing the field-vector from below.
            this.rotationSense = spin/Math.abs(spin) // == -1 if viewing from below

            /* PERIODICITY FROM SPIN?
             Why not simply use the accumulated spin-angle and time-span to calculate the rotation-period?
             --Because the arbitrary start and end angles may be subject to fore-shortening, and hence
             quite inaccurate! The only dependable angles occur as we pass over the major-axis (see below).
            */

            // The ratio of the extreme axis lengths gives the eccentricity of this Ellipse
            this.eccentricity = longest / shortest
            // Readings on a near-circular Ellipse won't ever be fore-shortened, so we can skip correction!
            this.isCircular = (this.eccentricity < Circular)

            /* We are trying to find a good approximation to the tilt of the Ellipse's major-axis.  
            We could simply nominate the longest candidate detected, but instead we will average them. 
            With an eccentric Ellipse, passing the major-axis twice per Spin-circle revolution yields neatly
            alternating candidates with "opposite" angles.
            With a more-nearly circular Ellipse, noisy readings can yield alternating clusters of local maxima 
            near each end of the axis.
            An almost circular Ellipse has no meaningful axis and generally yields multiple spurious candidates. 
            */

            // purge any local maximum whose vector length is too short --it's nowhere near the major-axis!
            let long = longest * LongEnough
            for (let i = 0; i < majors.length; i++) {
                if (majors[i].size < long) {
                    majors.splice(i, 1)  // disqualified!
                    i-- // (all subsequent candidates now shuffle up by one place!)
                }
            }
            // A simple way to form a consensus angle is to use vector addition
            let turns = 0 
            let endTime = 0
            let flipped = false
            let uSum = 0
            let vSum = 0
            let count = majors.length
            if (count > 0) {
                // the first candidate fixes which "end" of the axis we're choosing
                let front = majors[0].angle
                uSum = majors[0].u
                vSum = majors[0].v
                let startTime = majors[0].time
                for (let i = 1; i < count; i++) {
                    // get angle difference, as +/- 2pi
                    let deviate = ((ThreePi + majors[i].angle - front) % TwoPi) - Math.PI
                    // does it point nearer to the front? ...or to the back?
                    if (Math.abs(deviate) < HalfPi) {
                        // add this next arrow directly to the chain (no need to flip)
                        uSum += majors[i].u
                        vSum += majors[i].v
                        // the first unflipped Arrow after one or more flipped ones clocks a new revolution
                        if (flipped) {
                            flipped = false
                            turns++
                            endTime = majors[i].time
                        }
                    } else { // flip this arrow before chaining it, as it's pointing the "wrong" way
                        flipped = true
                        uSum -= majors[i].u
                        vSum -= majors[i].v
                    } 
                    
                    if ((mode == Mode.Trace) || (mode == Mode.Debug)) {
                            datalogger.log(
                                datalogger.createCV("plane", this.plane),
                                datalogger.createCV("time", majors[i].time),
                                datalogger.createCV("uAxis", majors[i].u),
                                datalogger.createCV("vAxis", majors[i].v),
                                datalogger.createCV("[angle]", asDegrees(majors[i].angle)),
                                datalogger.createCV("size", majors[i].size))
                    }
                }
                // normalise the resultant's vector coordinates (not strictly necessary!)
                uSum /= count
                vSum /= count
                // compute the average rotation time (so long as we've made at least one complete revolution)
                if (endTime > 0) { 
                    this.period = (endTime - startTime) / turns
                } else {
                    this.period = -1
                }
            } // else uSum & vSum remain at zero
            // create an Arrow pointing to {uSum,vSum}, showing the overall direction of the chained source Arrows
            this.majorAxis = new Arrow(uSum, vSum, 0)
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
    let north = 0 // reading registered as "North"
    let strength = 0 // the average magnetic field-strength observed by the magnetometer
    let period = -1 // overall assessment of average rotation time

    // correction parameters adopted from bestView Ellipse for future readings
    let rotationSense = 1 // set to -1 if orientation means field-vector projection is "from below"
    let isCircular: boolean // if bestView Ellipse is circular, no correction is needed
    let uOff: number // horizontal origin offset
    let vOff: number // vertical origin offset
    let theta: number  // major-axis tilt angle (in radians anticlockwise from the U-axis)
    let cosTheta: number; // saved for efficiency
    let sinTheta: number; //      ditto
    let scale: number // stretch-factor for correcting foreshortened readings (= eccentricity)

    // test-related globals
    let mode:Mode = Mode.Normal // mode switch for logging
    let dataset: string = "" // test dataset to use
    let testData: number[][] = [] //[X,Y,Z] magnetometer values for test cases
    let test = 0 // global selector for test-cases

    // EXPORTED USER INTERFACES   

    /** 
     * Assuming the buggy is currently spinning clockwise on the spot, capture a 
     * time-stamped sequence of magnetometer readings from which to set up the compass.
     *
     * @param ms scanning-time in millisecs (long enough for more than one full rotation)    
     */

    //% block="scan clockwise for (ms) $ms" 
    //% inlineInputMode=inline 
    //% ms.shadow="timePicker" 
    //% ms.defl=0 
    //% weight=90 

    export function scanClockwise(ms: number) {
        // Every 25-30 ms over the specified duration (generally a couple of seconds),
        // magnetometer readings are sampled and a new [X,Y,Z] triple added to the scanData[] array.
        // A timestamp for each sample is also recorded in the scanTimes[] array.

        // NOTE: to smooth out jitter, each reading is always a rolling sum of SEVEN consecutive
        // readings, effectively amplifying seven-fold the dynamic field range due to the Earth 
        // (from ~50 microTeslas to ~350)
        scanTimes = []
        scanData = []

        if (mode != Mode.Normal) {
            datalogger.deleteLog()
            datalogger.includeTimestamp(FlashLogTimeStampFormat.Milliseconds)
        }

        if (mode == Mode.Debug) {
            simulateScan(dataset)
            basic.pause(ms)
        } else { // use live magnetometer
            let index = 0
            let xRoll: number[] = []
            let yRoll: number[] = []
            let zRoll: number[] = []
            let x = 0
            let y = 0
            let z = 0
            let field = 0
            // discard the very first magnetometer readings as these are sometimes dubious!
            field = input.magneticForce(0)
            field = input.magneticForce(1)
            field = input.magneticForce(2)
            
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
            // while ((scanTimes.length < TooManySamples) && (now < finish)) {
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


        if ((mode == Mode.Trace) || (mode == Mode.Capture)) {
            datalogger.setColumnTitles("index", "t", "x", "y", "z")
            for (let i = 0; i < scanTimes.length; i++) {
                datalogger.log(
                    datalogger.createCV("index", i),
                    datalogger.createCV("t", scanTimes[i]),
                    datalogger.createCV("x", round2(scanData[i][Dimension.X])),
                    datalogger.createCV("y", round2(scanData[i][Dimension.Y])),
                    datalogger.createCV("z", round2(scanData[i][Dimension.Z])))
            }
        }
    }



    /**
     * Analyse the scanned data to prepare for reading compass-bearings.
     * Then read the magnetometer and register the buggy's current direction as "North",
     * (i.e. the direction that will in future return zero as its heading).
     * 
     * The actual direction of the buggy when this function is called is arbitrary:
     * it could be Magnetic North; or True North (compensating for local declination); 
     * or any convenient direction from which to measure subsequent heading angles.
     * 
     * @return zero if successful, or a negative error code:
     *
     *      -1 : NOT ENOUGH SCAN DATA

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
        let xlo = 9999999
        let ylo = 9999999
        let zlo = 9999999
        let xhi = -9999999
        let yhi = -9999999
        let zhi = -9999999
        for (let i = 0; i < nSamples; i++) {
            xhi = Math.max(xhi, scanData[i][Dimension.X])
            yhi = Math.max(yhi, scanData[i][Dimension.Y])
            zhi = Math.max(zhi, scanData[i][Dimension.Z])
            xlo = Math.min(xlo, scanData[i][Dimension.X])
            ylo = Math.min(ylo, scanData[i][Dimension.Y])
            zlo = Math.min(zlo, scanData[i][Dimension.Z])
        }

        // get RMS field-strength
        let xField = (xhi - xlo) / 2
        let yField = (yhi - ylo) / 2
        let zField = (zhi - zlo) / 2
        strength = Math.sqrt((xField * xField) + (yField * yField) + (zField * zField))

        // Bail out early if the scan didn't properly detect the Earth's magnetic field,
        // (perhaps due to magnetic shielding)
        if (strength < MarginalField) {
            return -2 // "FIELD STRENGTH TOO WEAK"
        }
        // The means of the extremes give the central offsets
        let xOff = (xhi + xlo) / 2
        let yOff = (yhi + ylo) / 2
        let zOff = (zhi + zlo) / 2

        // re-centre all of the scanData samples, so eliminating "hard-iron" magnetic effects
        for (let i = 0; i < nSamples; i++) {
            scanData[i][Dimension.X] -= xOff
            scanData[i][Dimension.Y] -= yOff
            scanData[i][Dimension.Z] -= zOff
        }

        // create three Ellipse instances for analysing each possible view
        views.push(new Ellipse("XY", Dimension.X, Dimension.Y, xOff, yOff))
        views.push(new Ellipse("YZ", Dimension.Y, Dimension.Z, yOff, zOff))
        views.push(new Ellipse("ZX", Dimension.Z, Dimension.X, zOff, xOff))

        // For each View, perform the analysis of eccentricity and Ellipse tilt-angle
        views[View.XY].analyseView()
        views[View.YZ].analyseView()
        views[View.ZX].analyseView()

        // check that at least one View saw at least one complete rotation (with a measurable period)...
        if ((views[View.XY].period == -1)
            && (views[View.YZ].period == -1)
            && (views[View.ZX].period == -1)) {
            period = -1
            return -3 // "NOT ENOUGH SCAN ROTATION"
        }

        // Choose the "roundest" Ellipse  --the one with lowest eccentricity.
        bestView = View.XY
        if (views[View.YZ].eccentricity < views[bestView].eccentricity) bestView = View.YZ
        if (views[View.ZX].eccentricity < views[bestView].eccentricity) bestView = View.ZX
/*        bestView = View.XY  // while debugging, force use of this view!
*/       
        // periodicity is unreliable in a near-circular View: average just the other two Views' measurements
        period = (views[0].period + views[1].period + views[2].period - views[bestView].period) / 2

        // For efficiency, extract various characteristics from the adopted "bestView" Ellipse
        uDim = views[bestView].uDim
        vDim = views[bestView].vDim
        uOff = views[bestView].uOff
        vOff = views[bestView].vOff
        scale = views[bestView].eccentricity
        theta = views[bestView].majorAxis.angle // the rotation (in radians) of the major-axis from the U-axis
        cosTheta = Math.cos(theta)
        sinTheta = Math.sin(theta)
        isCircular = views[bestView].isCircular
        rotationSense = views[bestView].rotationSense

        // Having successfully set up the projection parameters for the bestView, get a
        // stable fix on the current heading, which we will then designate as "North".
        // (This is the global fixed bias to be subtracted from all future readings)
        north = takeSingleReading()

        if ((mode == Mode.Trace) || (mode == Mode.Debug)) {
            datalogger.log(
                datalogger.createCV("view", views[bestView].plane),
                datalogger.createCV("scale", scale),
                datalogger.createCV("circular?", isCircular),
                datalogger.createCV("sense?", rotationSense),
                datalogger.createCV("theta", theta),
                datalogger.createCV("[theta]", round2(asDegrees(theta))),
                datalogger.createCV("north", north),
                datalogger.createCV("[north]", round2(asDegrees(north))),
                datalogger.createCV("period", period),
            )
        }

        // we've now finished with the scanning data and Ellipse objects, so release their memory
        scanTimes = []
        scanData = []
        views = []

        // SUCCESS!
        return 0
    }


    /**
     * Read the magnetometer
     * 
     * @return the current heading of the buggy
     * 
     * (in degrees clockwise, relative to "North")
     */
    //% block="degrees" 
    //% inlineInputMode=inline 
    //% weight=70
    export function degrees(): number {
    // Depending on mounting orientation, the bestView might possibly be seeing the
    // Spin-Circle from "underneath", with the field-vector appearing to move clockwise 
    // --effectively experiencing an anti-clockwise scan. 
    // In this case the rotationSense will be negative.
        return asDegrees((takeSingleReading() - north) * rotationSense)
        // NOTE: that there is a double reversal going on here:
        // Viewed from above, the Field-vector reading increases (anticlockwise) w.r.t "North"
        // as the buggy's heading increases (clockwise). From below, a third reversal is needed!
    }

    /**
     * The average rotation time of the most recent scan 
     * @return rotation time, or error-value:
     * 
     *      -4 : SUCCESSFUL SCAN IS NEEDED FIRST
     */
    //% block="spin time (ms)" 
    //% inlineInputMode=inline 
    //% weight=60 
    export function spinTime(): number {
        if (period == -1) {
            return -4 // ERROR: SUCCESSFUL SCAN IS NEEDED FIRST
        } else {
            return period
        }
    }

    /**
     * The average rotation rate of the most recent scan 
     * 
     * @return revs-per-minute, or error value:
     * 
     *      -4 : SUCCESSFUL SCAN IS NEEDED FIRST
     */
    //% block="spin rate (RPM)" 
    //% inlineInputMode=inline 
    //% weight=50 
    export function spinRate(): number {
        if (period == -1) {
            return -4 // ERROR: SUCCESSFUL SCAN IS NEEDED FIRST
        } else {
            return 60000 / period
        }
    }

    /**
     * While scanning, wheels are rotated in opposite directions, giving a spin-rate for the 
     * selected power setting. Based on the axle-length and latest spin-rate, this function 
     * estimates the forward speed to be expected when using that power setting.
     * (NOTE that tyre-friction or skidding when turning may make this a fairly inaccurate estimate!)
     * 
     * @param axleLength : distance betweeen mid-lines of tyres (in mm)
     * 
     * @return speed in mm-per-second, or error value:
     * 
     *      -4 : SUCCESSFUL SCAN IS NEEDED FIRST
     */
    //% block="equivalent speed (mm/s), axle length (mm) = $axleLength" 
    //% inlineInputMode=inline 
    //% weight=50 
    export function equivalentSpeed(axleLength: number): number {
        if ((views.length == 0) || (views[bestView].period <= 0)) {
            return -4 // ERROR: SUCCESSFUL SCAN IS NEEDED FIRST
        } else {
            // compute tangential speed of wheel-centre in mm/s:
            // it takes [period] ms to cover [2pi * axleLength/2] mm
            return (Math.PI * axleLength * 1000 / period)
        }
    }

  

    /**
     * Choose mode: whether to run normally,
     *  - or to use Data Logger to grab a new test dataset,
     *  - or to debug processing using a named test dataset
     * (new datasets need editing externally before compiling into simulateScan)
     */
    //% block="reset to mode $newMode using $name"
    //% inlineInputMode=inline 
    //% weight=10
    export function resetMode(newMode:Mode, name: string) {
        // always reinitialise key data
        scanTimes = []
        scanData = []
        views = []
        bestView = -1
        uDim = -1
        vDim = -1
        period = -1
        north = 0
        testData = []
        test = 0
        // now adopt the new mode
        mode = newMode
        dataset = name
    }

    // UTILITY FUNCTIONS

    /** Take the sum of seven new readings to get a stable fix on the current heading.
     *  @return projected angle of the magnetic field-vector (in radians anticlockwise
     * from the horizontal U-axis)
     */

    /* Although eventually we'd only need [uDim, vDim], we'll sum and log all three Dims.
       This will allow us, while testing, to override automatic choice of bestView
       and check out more severe levels of correction! Eventually, optimise this out!
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
        if (mode == Mode.Debug) { // just choose the next test-data value
            uRaw = testData[test][uDim]
            vRaw = testData[test][vDim]
            test = (test + 1) % testData.length
        } else {
            let xyz = [0, 0, 0]

            // get a new sample as the sum of seven readings, 10ms apart
            xyz[0] = input.magneticForce(0)
            xyz[1] = input.magneticForce(1)
            xyz[2] = input.magneticForce(2)

            for (let i = 0; i < 6; i++) {
                basic.pause(10)
                xyz[0] += input.magneticForce(0)
                xyz[1] += input.magneticForce(1)
                xyz[2] += input.magneticForce(2)
            }

            if ((mode == Mode.Trace)||(mode == Mode.Capture)) {
                datalogger.log(
                    datalogger.createCV("index", test),
                    datalogger.createCV("x", round2(xyz[0])),
                    datalogger.createCV("y", round2(xyz[1])),
                    datalogger.createCV("z", round2(xyz[2])))
            }

            // even in normal operation, it's currently useful to keep a history of readings
            testData.push(xyz)
            test++ // clock another test sample

            // use just our current bestView's two dimensions
            uRaw = xyz[uDim]
            vRaw = xyz[vDim]
        }

        // re-centre this latest point w.r.t our Ellipse origin
        u = uRaw - uOff
        v = vRaw - vOff

        if (isCircular) {
            reading = Math.atan2(v, u) 
        } else { 
        // Unless this Ellipse.isCircular, any {u,v} reading will be foreshortened in this View, and
        // must be stretched along the Ellipse minor-axis to place it correctly onto the Spin-Circle.

            // First rotate CLOCKWISE by theta (so aligning the Ellipse minor-axis angle with the V-axis)
            uNew = u * cosTheta + v * sinTheta
            vNew = v * cosTheta - u * sinTheta
            // Now scale up along V, re-balancing the axes to make the Ellipse circular
            uFix = uNew
            vFix = vNew * scale
            // get the adjusted angle for this corrected {u,v}
            reading = Math.atan2(vFix, uFix)
            // finally, undo the rotation by theta
            reading += theta
        }
        
        if ((mode == Mode.Trace)||(mode == Mode.Debug)) {
            // just for debug, show coordinates of "stretched" reading after undoing rotation
            let uStretch = uFix * cosTheta - vFix * sinTheta
            let vStretch = vFix * cosTheta + uFix * sinTheta
            datalogger.log(
                datalogger.createCV("u", round2(u)),
                datalogger.createCV("v", round2(v)),
                datalogger.createCV("uNew", round2(uNew)),
                datalogger.createCV("vNew", round2(vNew)),
                datalogger.createCV("uFix", round2(uFix)),
                datalogger.createCV("vFix", round2(vFix)),
                datalogger.createCV("uStretch", round2(uStretch)),
                datalogger.createCV("vStretch", round2(vStretch)),
                datalogger.createCV("reading", round2(reading)),
                datalogger.createCV("[reading]", round2(asDegrees(reading) * rotationSense))
            )
        }
        return reading
    }

    // helpful for logging...
    function round2(v: number): number {
        return (Math.round(100 * v) / 100)
    }


    // Convert an angle measured in radians to degrees.
    function asDegrees(angle: number): number {
        return ((angle * RadianDegrees) + 360) % 360
    }

    // While debugging, it is necessary to re-use predictable sample data for a variety of use-cases 
    // captured (using the datalogger) from earlier live runs. 
    // [This function is greedy on memory and can be commented-out once the extension has been fully debugged.]
    function simulateScan(dataset: string) {
        let xData: number[] = []
        let yData: number[] = []
        let zData: number[] = []
        let xTest: number[] = []
        let yTest: number[] = []
        let zTest: number[] = []
        switch (dataset) {
/*
            case "angled": // mounted at strange angle with approx spin axis components [2,-2,1]; dip=70
                scanTimes = [15845, 15865, 15885, 15905, 15925, 15945, 15965, 15985, 16005, 16025, 16045, 16065, 16085, 16105, 16125, 16145, 16165, 16185, 16207, 16229, 16249, 16269, 16289, 16309, 16329, 16349, 16369, 16389, 16409, 16429, 16449, 16469, 16489, 16509, 16529, 16549, 16569, 16589, 16609, 16629, 16649, 16669, 16689, 16709, 16729, 16749, 16769, 16789, 16809, 16829, 16849, 16869, 16889, 16909, 16929, 16949, 16969, 16989, 17009, 17029, 17049, 17069, 17089, 17109, 17129, 17153, 17173, 17193, 17213, 17233, 17253, 17273, 17293, 17313, 17333, 17353, 17373, 17393, 17413, 17433, 17453, 17473, 17493, 17513, 17533, 17553, 17573, 17593, 17613, 17633, 17653, 17673, 17693, 17713, 17733, 17753, 17773, 17793, 17813, 17837, 17857, 17877, 17897, 17917, 17937, 17957, 17977, 17997, 18017, 18037, 18057, 18077, 18097, 18117, 18137, 18157, 18177, 18197, 18221, 18241, 18261, 18281, 18301, 18321, 18341, 18361, 18381, 18401, 18421, 18441, 18461, 18481, 18501, 18521, 18541, 18561, 18581, 18601, 18621, 18644, 18665, 18685, 18705, 18725, 18745, 18765, 18785, 18805, 18825, 18845, 18865, 18885, 18905, 18925, 18945, 18965, 18989, 19009, 19029, 19049, 19069, 19089, 19109, 19129, 19149, 19169, 19189, 19209, 19233, 19253, 19273, 19293, 19313, 19333, 19353, 19373, 19397, 19421, 19441, 19461, 19481, 19501, 19521, 19541, 19561, 19581, 19601, 19621, 19641, 19661, 19681, 19701, 19721, 19741, 19761, 19781, 19801, 19821, 19846, 19869, 19889, 19909, 19929, 19949, 19969, 19989, 20009, 20029, 20049, 20069, 20089, 20109, 20129, 20149, 20169, 20189, 20209, 20229, 20249, 20269, 20289, 20309, 20329, 20349, 20369, 20389, 20409, 20429, 20449, 20469, 20489, 20509, 20529, 20549, 20569, 20589, 20609, 20629, 20649, 20669, 20689, 20709, 20729, 20749, 20769, 20789, 20809, 20829, 20849, 20869, 20889, 20909, 20929, 20949, 20969, 20997, 21017, 21037, 21057, 21077, 21097, 21117, 21137, 21157, 21177, 21197, 21217, 21237, 21257, 21277, 21297, 21317, 21337, 21357, 21377, 21397, 21417, 21437, 21457, 21477, 21497, 21517, 21537, 21557, 21577, 21597, 21617, 21637, 21657, 21677, 21697, 21717, 21737, 21757, 21785, 21805, 21825, 21845, 21865, 21885, 21905, 21925, 21945, 21965, 21985, 22005, 22025, 22045, 22065, 22085, 22105, 22125, 22145, 22165, 22185, 22205, 22225, 22245, 22265, 22285, 22305, 22325, 22353, 22373, 22393, 22413, 22433, 22453, 22473, 22493, 22513, 22533, 22553, 22573, 22593, 22613, 22633, 22660, 22681, 22701, 22721]
                xData = [-253.65, -254.4, -254.1, -253.5, -253.2, -253.35, -252.9, -251.25, -249.45, -248.4, -247.05, -244.8, -241.95, -237.75, -234.75, -230.85, -226.95, -223.2, -218.85, -214.05, -210.15, -205.05, -199.5, -193.2, -185.85, -178.35, -172.65, -165.75, -158.1, -152.4, -146.1, -140.7, -135.45, -128.4, -122.25, -117.3, -111.15, -105.9, -100.8, -96.3, -91.65, -87.45, -83.25, -78.75, -73.2, -69.45, -65.25, -62.25, -59.25, -56.25, -53.85, -52.2, -49.5, -46.65, -44.55, -42, -39.3, -37.95, -36.45, -35.1, -34.65, -34.05, -33.6, -33.45, -32.55, -32.85, -32.85, -32.7, -31.5, -32.55, -34.05, -34.65, -35.1, -36.3, -37.8, -39.9, -40.5, -42.6, -45.75, -47.55, -49.2, -50.85, -54, -55.8, -56.85, -57.75, -61.05, -64.5, -67.95, -70.35, -74.55, -78.3, -82.65, -85.5, -89.1, -92.55, -95.7, -98.7, -101.4, -103.8, -108.3, -111.45, -114.9, -119.1, -122.85, -127.05, -131.7, -134.85, -139.05, -142.95, -145.8, -148.5, -151.2, -154.35, -157.8, -160.5, -164.4, -169.35, -174.3, -179.4, -183.45, -186.45, -190.2, -193.35, -196.8, -199.8, -202.5, -205.95, -209.7, -213.15, -214.8, -217.5, -220.5, -223.05, -224.4, -226.5, -228.6, -232.65, -234.9, -238.05, -241.05, -244.8, -246.75, -247.8, -249.15, -250.5, -250.5, -251.7, -252.3, -254.4, -256.35, -258, -259.05, -259.95, -260.85, -260.4, -259.05, -258.15, -256.2, -254.7, -253.5, -251.4, -250.8, -249.6, -247.95, -246.3, -243.9, -241.35, -238.05, -234.3, -230.4, -227.1, -224.25, -219.45, -215.25, -211.8, -208.2, -204.6, -199.35, -193.95, -190.65, -185.55, -179.1, -171.3, -164.25, -159, -152.7, -145.5, -140.1, -134.25, -130.35, -125.25, -119.25, -113.4, -108.45, -103.65, -99.3, -94.65, -91.35, -87.3, -84.6, -81.45, -77.25, -73.5, -69.6, -65.7, -63.15, -58.8, -55.5, -52.65, -50.85, -49.05, -46.95, -45.6, -44.25, -42.3, -40.5, -38.7, -36.45, -35.7, -33.9, -33.15, -33.3, -33.6, -33.6, -34.65, -34.95, -35.7, -36.3, -36.75, -37.35, -37.95, -38.55, -38.7, -38.85, -39.75, -39.75, -40.65, -41.55, -43.8, -45.45, -47.55, -49.5, -51.45, -52.95, -54.45, -55.5, -57.75, -58.5, -60.6, -63.3, -65.85, -68.25, -70.2, -72.45, -75.9, -77.55, -80.1, -82.2, -84.6, -87.6, -89.1, -91.95, -95.1, -97.65, -100.5, -104.25, -107.55, -110.7, -113.1, -115.65, -118.5, -121.35, -124.05, -125.55, -128.55, -130.65, -132.9, -135.6, -137.55, -139.35, -142.8, -145.5, -148.05, -150.45, -152.4, -154.95, -157.5, -159.15, -161.25, -163.65, -166.65, -169.65, -172.5, -175.5, -178.95, -182.55, -185.7, -188.7, -191.55, -194.4, -198.15, -200.85, -203.85, -207.6, -210.3, -212.7, -216, -217.2, -219.45, -221.55, -222.9, -225.15, -227.25, -228.75, -230.25, -232.65, -234.15, -236.7, -238.35, -240.3, -241.8, -244.8, -245.7, -246.75, -246.15, -247.65, -248.7, -250.5, -250.8, -250.5, -251.4, -252.6, -253.5, -254.1, -253.65, -253.05, -253.95, -253.35, -253.65, -252]
                yData = [300.75, 304.2, 310.05, 314.1, 318.6, 323.25, 327.6, 331.35, 336.6, 339.75, 343.95, 349.05, 354.3, 357.75, 363.3, 366.3, 369.45, 374.1, 376.5, 378.3, 381, 382.5, 385.5, 387.75, 388.65, 389.85, 390.75, 391.35, 392.4, 391.65, 390.75, 390.15, 388.95, 386.85, 385.05, 383.1, 381.3, 379.35, 376.65, 374.7, 372.6, 370.05, 366.6, 363.45, 360.15, 355.8, 351.75, 348.9, 344.85, 340.2, 337.5, 333.6, 331.05, 327.45, 323.4, 320.25, 316.95, 312.15, 307.8, 302.25, 297.3, 292.35, 286.5, 282.6, 277.05, 273.15, 269.55, 266.1, 261.6, 256.95, 251.7, 247.65, 243.75, 239.4, 234.45, 231.15, 228.45, 225.9, 222.9, 219.45, 216.75, 214.65, 211.65, 209.25, 206.4, 204.15, 202.5, 200.55, 198.75, 196.5, 195, 193.35, 191.55, 190.05, 188.55, 186.6, 185.7, 185.25, 184.35, 183.75, 182.85, 181.95, 181.95, 182.1, 181.8, 182.25, 182.4, 182.55, 183.9, 185.1, 185.1, 184.65, 186.15, 187.65, 188.7, 189.45, 191.55, 193.8, 195, 196.05, 197.25, 198.75, 200.7, 201.3, 203.4, 206.7, 209.1, 211.2, 214.8, 217.5, 221.1, 223.35, 224.7, 227.25, 230.25, 231.9, 234.9, 236.7, 238.95, 243.15, 246.6, 250.2, 253.65, 257.1, 261.6, 266.1, 271.05, 274.35, 278.4, 282.9, 286.95, 290.25, 294.3, 297.15, 300.6, 304.2, 307.95, 311.25, 314.7, 318.45, 321.9, 326.4, 330.3, 333.9, 338.1, 343.05, 346.65, 350.4, 353.4, 357.15, 361.65, 365.85, 368.25, 372, 375.6, 379.35, 381.9, 384, 386.55, 388.95, 390.15, 391.8, 392.4, 393.6, 394.05, 392.85, 391.95, 391.2, 390.3, 389.1, 387.45, 385.35, 383.85, 382.65, 381, 378.45, 375.6, 371.7, 369.9, 367.95, 364.95, 361.65, 358.95, 356.7, 354.6, 350.25, 345.45, 341.55, 337.95, 334.35, 330.15, 326.25, 322.65, 318.45, 315, 311.4, 307.5, 303.75, 300, 295.2, 291.6, 286.5, 282.3, 277.8, 273.9, 270.15, 267.75, 264.75, 261.75, 259.2, 255.75, 253.95, 250.5, 245.85, 241.8, 238.8, 234.9, 232.2, 228.3, 225.6, 223.65, 221.55, 218.85, 217.35, 215.1, 212.25, 210, 208.8, 207.15, 205.95, 203.7, 202.5, 200.4, 199.5, 197.7, 196.05, 194.85, 193.35, 191.55, 191.4, 190.2, 189.3, 187.65, 186, 184.35, 183.45, 181.8, 181.35, 180.15, 181.05, 181.8, 182.7, 183.15, 184.5, 184.35, 184.05, 184.35, 183.6, 184.05, 184.65, 184.5, 185.1, 186.9, 186.6, 186.75, 186.6, 186, 186.15, 187.2, 187.05, 188.7, 190.05, 190.8, 192.6, 193.95, 194.25, 196.2, 197.4, 199.8, 201.15, 202.8, 204.75, 206.7, 208.05, 209.25, 211.05, 213.75, 216.3, 218.55, 221.1, 223.8, 226.35, 228.6, 231.15, 233.1, 235.65, 237.9, 240.3, 242.85, 245.1, 248.55, 252.45, 255.45, 258.6, 261.45, 265.35, 268.95, 270.9, 273, 276.6, 280.5, 284.7, 288.3, 292.2, 297.3, 301.35, 304.5, 307.05, 309.9, 312.75, 316.65]
                zData = [442.2, 439.65, 438.15, 436.2, 432.6, 429.75, 427.65, 424.5, 421.5, 418.35, 415.95, 414, 411.6, 408.45, 405.75, 402.9, 398.85, 395.85, 391.8, 387.9, 383.7, 379.65, 376.05, 371.85, 368.1, 365.4, 361.65, 358.8, 355.8, 353.1, 351.6, 349.2, 345.75, 344.25, 342.3, 341.55, 340.35, 339.3, 338.25, 338.55, 338.25, 337.8, 336.9, 336.45, 336.6, 335.7, 334.8, 335.1, 335.4, 336, 335.7, 336.9, 338.4, 339.6, 340.35, 342.15, 343.95, 345.6, 346.65, 349.65, 351.6, 353.7, 355.35, 356.85, 358.8, 361.2, 362.1, 364.8, 366.45, 369.6, 370.8, 373.65, 375.6, 377.1, 379.05, 382.5, 384.3, 388.5, 390.6, 393.3, 396.6, 398.55, 400.2, 402.75, 403.95, 406.65, 408.15, 410.55, 412.5, 414.3, 416.1, 419.4, 421.35, 423.75, 425.25, 427.35, 429.75, 432.6, 434.1, 435.9, 438.3, 439.8, 441.3, 443.25, 444.15, 444.75, 446.4, 447.3, 448.8, 450.75, 451.65, 451.95, 453.75, 454.05, 454.35, 455.4, 456, 456.75, 457.2, 456.9, 457.65, 459.3, 460.2, 461.1, 461.4, 462.15, 462.75, 462.6, 462, 462.6, 462.75, 462.75, 463.35, 463.8, 465.3, 465.3, 464.1, 463.35, 462.75, 461.4, 460.5, 457.8, 456, 455.25, 453.45, 451.2, 450, 448.8, 447.3, 446.4, 444.15, 442.8, 441.75, 440.4, 437.4, 436.5, 435.15, 433.5, 432.15, 431.55, 429.3, 427.5, 425.25, 421.8, 419.1, 417, 413.4, 409.95, 407.1, 403.2, 400.8, 397.65, 392.7, 388.5, 386.1, 382.8, 379.95, 377.4, 374.55, 372.3, 369.45, 367.2, 363.75, 361.35, 358.05, 354.75, 351.3, 349.5, 346.5, 345.45, 342.6, 340.5, 339.3, 339.9, 338.7, 337.2, 336.15, 336.6, 336, 336.3, 335.1, 335.4, 336.15, 336.3, 336.45, 337.2, 337.8, 338.25, 339.15, 339.3, 340.65, 341.55, 341.4, 342.45, 343.05, 343.65, 345.9, 346.95, 348.6, 350.7, 351.3, 353.4, 355.35, 356.85, 357.9, 358.8, 360.45, 363.45, 364.95, 366.6, 368.85, 372.15, 374.55, 376.8, 377.25, 380.1, 381.45, 382.5, 382.95, 384.9, 387.9, 391.5, 393.15, 394.95, 396.6, 399.9, 400.95, 402, 402.75, 404.25, 406.5, 408.6, 409.35, 411.45, 412.95, 414.3, 416.4, 418.2, 418.8, 420.75, 421.35, 423.9, 424.95, 426.15, 427.95, 430.2, 432, 435, 434.85, 438.15, 439.8, 439.65, 441.75, 441.6, 442.35, 444.9, 444.9, 445.65, 448.2, 448.8, 449.55, 450, 450.3, 450.75, 453.15, 453.6, 453.45, 454.8, 455.55, 456.15, 456.9, 455.4, 456.15, 457.05, 458.1, 458.85, 459.75, 460.8, 462.9, 462.6, 463.5, 464.1, 464.55, 464.4, 463.5, 462.45, 463.05, 462.75, 462.45, 462.45, 461.7, 462.9, 463.05, 462, 462, 462.45, 461.1, 460.5, 460.35, 459.15, 459.9, 459, 457.2, 456.75, 456, 454.2, 454.05, 452.7, 451.05, 449.4, 448.05, 447, 446.25, 443.55, 441.75, 440.25, 440.25, 438.75, 436.65, 433.35, 432.15]
                xTest = [-49.05, -47.25, -34.95, -79.8, -163.8, -231.15, -252.75, -201.6, -120.15, -49.05]
                yTest = [333.9, 330.3, 253.95, 191.4, 191.85, 242.55, 324.15, 385.5, 383.7, 328.95]
                zTest = [336.75, 338.1, 368.1, 412.8, 453.9, 453.75, 423, 368.85, 336.3, 333.3]
                break

            case "zdown70": // mounted face-up so Z-axis pointing down; dip=70
                scanTimes = [48288, 48308, 48328, 48348, 48368, 48388, 48408, 48428, 48448, 48468, 48488, 48508, 48528, 48548, 48568, 48588, 48608, 48628, 48650, 48672, 48692, 48712, 48732, 48752, 48772, 48792, 48812, 48832, 48852, 48872, 48892, 48912, 48932, 48952, 48972, 48992, 49012, 49032, 49052, 49072, 49092, 49112, 49132, 49152, 49172, 49192, 49212, 49232, 49252, 49272, 49292, 49312, 49332, 49352, 49372, 49392, 49412, 49432, 49452, 49472, 49492, 49512, 49532, 49552, 49572, 49596, 49616, 49636, 49656, 49676, 49696, 49716, 49736, 49756, 49776, 49796, 49816, 49836, 49856, 49876, 49896, 49916, 49936, 49956, 49976, 49996, 50016, 50036, 50056, 50076, 50096, 50116, 50136, 50156, 50176, 50196, 50216, 50236, 50256, 50280, 50300, 50320, 50340, 50360, 50380, 50400, 50420, 50440, 50460, 50480, 50500, 50520, 50540, 50560, 50580, 50600, 50620, 50640, 50663, 50684, 50704, 50724, 50744, 50764, 50784, 50804, 50824, 50844, 50864, 50884, 50904, 50924, 50944, 50964, 50984, 51004, 51024, 51044, 51064, 51084, 51108, 51128, 51148, 51168, 51188, 51208, 51228, 51248, 51268, 51288, 51308, 51328, 51348, 51368, 51388, 51408, 51431, 51452, 51472, 51492, 51512, 51532, 51552, 51572, 51592, 51612, 51632, 51652, 51675, 51696, 51716, 51736, 51756, 51776, 51796, 51816, 51836, 51860, 51880, 51900, 51920, 51940, 51960, 51980, 52000, 52020, 52040, 52060, 52080, 52100, 52120, 52140, 52160, 52180, 52200, 52220, 52240, 52260, 52285, 52308, 52328, 52348, 52368, 52388, 52408, 52428, 52448, 52468, 52488, 52508, 52528, 52548, 52568, 52588, 52608, 52628, 52648, 52668, 52688, 52708, 52728, 52748, 52768, 52788, 52808, 52828, 52848, 52868, 52888, 52908, 52928, 52948, 52968, 52988, 53008, 53028, 53048, 53068, 53088, 53108, 53128, 53148, 53168, 53188, 53208, 53228, 53248, 53268, 53288, 53308, 53328, 53348, 53368, 53388, 53408, 53428, 53454, 53476, 53496, 53516, 53536, 53556, 53576, 53596, 53616, 53636, 53656, 53676, 53696, 53716, 53736, 53756, 53776, 53796, 53816, 53836, 53856, 53876, 53896, 53916, 53936, 53956, 53976, 53996, 54016, 54036, 54056, 54076, 54096, 54116, 54136, 54156, 54176, 54196, 54224, 54244, 54264, 54284, 54304, 54324, 54344, 54364, 54384, 54404, 54424, 54444, 54464, 54484, 54504, 54524, 54544, 54564, 54584, 54604, 54624, 54644, 54664, 54684, 54704, 54724, 54744, 54771, 54792, 54812, 54832, 54852, 54872, 54892, 54912, 54932, 54952, 54972, 54992, 55012, 55032, 55052, 55072, 55099, 55120, 55140, 55160, 55180]
                xData = [-183.45, -191.1, -194.25, -197.1, -199.95, -203.1, -206.25, -209.25, -212.7, -216.9, -220.95, -224.4, -228.45, -231.6, -235.5, -239.1, -243.3, -246.15, -250.5, -253.65, -256.8, -260.1, -263.4, -266.85, -270.45, -273.45, -275.85, -279, -283.2, -286.95, -289.65, -293.7, -297.75, -303.3, -307.65, -310.95, -315.3, -320.7, -325.8, -330.45, -335.4, -340.5, -344.25, -347.7, -351.9, -354.75, -357.75, -359.4, -362.1, -366.45, -370.05, -371.55, -373.35, -374.85, -376.8, -377.7, -377.25, -377.1, -377.7, -378.15, -378, -376.95, -375.15, -375.15, -373.65, -370.05, -367.2, -363.9, -359.55, -355.8, -350.1, -344.7, -340.05, -334.5, -329.55, -325.5, -320.25, -314.85, -308.55, -303.9, -298.2, -292.2, -285.45, -279.9, -274.05, -268.95, -262.65, -256.8, -250.95, -246, -240.75, -235.2, -230.4, -225, -219, -214.35, -208.95, -204.15, -199.8, -194.85, -189.3, -184.35, -178.95, -173.4, -168.15, -163.2, -158.4, -155.25, -152.4, -148.5, -146.7, -144, -141.75, -140.7, -138.6, -138.3, -139.05, -138.6, -139.2, -139.65, -139.2, -139.5, -140.25, -140.4, -141.45, -142.5, -143.25, -144.6, -147.3, -148.8, -151.35, -154.2, -156.9, -160.95, -165.15, -169.5, -173.1, -177.6, -180.3, -184.65, -188.25, -192, -195.75, -200.4, -204.75, -210.15, -214.05, -217.95, -222.75, -226.5, -231, -234.45, -238.65, -243, -247.95, -252.15, -257.1, -261.9, -267.3, -273, -278.55, -284.55, -290.25, -295.35, -299.4, -304.2, -308.1, -312.15, -316.35, -319.2, -323.1, -326.4, -329.1, -332.1, -334.8, -336.45, -340.65, -343.95, -348, -351.9, -355.5, -359.85, -364.5, -367.8, -370.5, -372.9, -375.45, -377.25, -378.15, -378.3, -378.6, -378.3, -379.35, -378.75, -378.6, -377.4, -376.5, -374.55, -373.2, -370.05, -367.65, -364.05, -360.9, -357.6, -354.45, -350.7, -346.2, -341.4, -337.35, -333.45, -328.8, -323.85, -319.65, -315.3, -310.35, -305.1, -300.3, -295.95, -290.55, -283.95, -278.85, -274.5, -270, -264.75, -258.45, -254.7, -250.35, -246.3, -240.6, -236.4, -231.6, -227.55, -222.3, -218.55, -213.3, -210.3, -204.6, -199.5, -194.85, -190.35, -184.8, -181.05, -175.95, -171.75, -168.15, -164.7, -159.75, -156.9, -153.45, -151.2, -149.55, -148.5, -146.25, -145.8, -144.3, -142.05, -140.1, -139.35, -137.85, -138.15, -137.7, -137.85, -138.45, -139.95, -139.8, -140.25, -140.7, -142.65, -144, -146.1, -146.4, -149.55, -151.8, -154.65, -156.3, -160.05, -163.8, -168.9, -172.05, -175.2, -179.1, -183.75, -187.05, -190.2, -192.9, -196.35, -201.45, -205.2, -209.25, -213.45, -217.5, -222.9, -228.3, -233.25, -238.65, -244.5, -250.35, -256.2, -261.15, -265.95, -271.2, -276.15, -280.35, -283.95, -288.9, -294.75, -299.85, -304.65, -308.85, -314.25, -319.5, -324.45, -326.7, -330.9, -334.8, -338.4, -341.55, -345.15, -348, -352.35, -355.05, -357.6, -361.35, -364.05, -367.05, -369.45, -371.25, -372.9, -374.85, -375.15, -376.8, -377.7, -377.85, -377.7, -378.15, -378.45, -378.6, -377.1, -375.9, -374.7, -373.05, -370.35]
                yData = [80.4, 81.3, 79.95, 76.5, 75.75, 74.7, 72.6, 70.65, 69.75, 68.25, 67.95, 65.4, 64.35, 63.75, 63.6, 63.6, 63.15, 62.1, 61.8, 61.8, 61.5, 60.3, 59.7, 59.1, 59.4, 60.9, 60.45, 61.5, 62.7, 62.7, 63.45, 64.2, 64.65, 66, 67.35, 69.75, 71.85, 74.55, 78, 82.05, 86.55, 89.7, 93.15, 97.5, 102.45, 106.95, 110.55, 115.35, 119.85, 124.35, 130.35, 135.15, 139.8, 145.8, 149.4, 155.25, 161.25, 166.8, 172.65, 179.4, 185.25, 192.3, 198.9, 204.75, 211.35, 217.2, 222.75, 228.15, 234.3, 240, 246.75, 251.1, 255.45, 260.4, 264.9, 268.5, 271.8, 273.3, 276, 279.45, 282, 284.1, 286.05, 286.95, 288.45, 289.5, 289.65, 289.65, 289.5, 288.9, 289.5, 289.5, 288.9, 287.7, 285.75, 284.4, 282.75, 279.15, 277.05, 273.3, 270.3, 266.1, 261.3, 256.5, 251.85, 246.75, 241.5, 235.8, 231.15, 226.95, 222, 217.65, 212.1, 207.6, 202.8, 197.7, 192, 186.9, 181.05, 175.95, 171, 166.65, 161.85, 157.8, 153, 148.95, 144.3, 139.95, 134.55, 130.05, 124.5, 119.55, 114, 108.9, 103.65, 99.6, 95.1, 91.2, 88.2, 85.5, 83.1, 81, 77.7, 74.85, 72, 70.2, 68.55, 66.6, 64.65, 63.45, 63.15, 63.15, 61.5, 60.75, 59.7, 59.7, 60.6, 60.3, 60, 60.9, 62.1, 63.75, 64.95, 64.8, 66.45, 68.25, 70.2, 70.8, 71.85, 73.2, 75.75, 77.1, 79.65, 81.6, 84.9, 88.05, 91.05, 93.9, 98.1, 102.3, 106.5, 111.45, 116.55, 120.75, 126.75, 132, 136.35, 142.35, 147.3, 153, 159.6, 165.45, 171.15, 177.75, 183.3, 189.3, 195.15, 201.15, 206.25, 211.95, 216.9, 222.75, 228, 233.25, 238.65, 243.45, 248.55, 253.95, 258.3, 262.8, 265.95, 269.4, 273.75, 276.45, 279, 280.95, 282.45, 286.5, 287.4, 288.45, 289.65, 289.95, 289.95, 291.45, 290.55, 291.75, 291.9, 291.6, 291.6, 292.2, 290.25, 289.65, 288.3, 286.65, 286.2, 283.8, 281.4, 279.9, 277.2, 274.05, 270.75, 266.85, 263.55, 259.5, 254.55, 250.05, 246.15, 242.25, 237.9, 234.15, 230.4, 226.5, 222, 217.35, 212.1, 207.15, 201.45, 195.9, 191.7, 185.55, 179.7, 174.75, 168.9, 163.2, 157.2, 151.05, 147, 141.15, 135.45, 130.2, 126.75, 122.25, 118.05, 114.3, 111.15, 107.85, 104.7, 100.8, 97.8, 94.35, 90.15, 87.45, 84.45, 82.35, 79.95, 78, 75.75, 74.4, 72.9, 70.95, 69.15, 66.9, 65.7, 64.95, 64.2, 62.55, 61.95, 61.35, 61.35, 61.35, 61.5, 61.5, 62.85, 64.8, 65.7, 67.05, 68.7, 70.05, 71.25, 73.05, 74.4, 76.95, 79.8, 82.8, 86.7, 90, 93.75, 96.9, 100.35, 103.65, 106.05, 109.35, 114.75, 119.4, 124.95, 130.05, 136.35, 142.35, 147.3, 151.8, 156.3, 161.4, 166.35, 170.7, 175.8, 180.9, 185.55, 190.5, 195.6, 201, 206.7, 211.5]
                zData = [420.6, 424.65, 424.65, 424.35, 424.65, 425.55, 425.25, 424.8, 424.2, 423.75, 423.6, 424.2, 424.8, 424.5, 424.05, 423.9, 424.05, 424.65, 422.85, 421.35, 420.9, 421.35, 420.6, 420.45, 418.35, 418.5, 418.8, 418.65, 418.2, 417.75, 418.35, 420.75, 420.75, 421.05, 421.2, 420.9, 421.05, 420.3, 418.8, 418.8, 417.6, 417.6, 417.6, 418.05, 417.6, 417.9, 417.75, 418.65, 417.9, 417, 416.7, 417.45, 416.55, 415.95, 415.2, 415.65, 415.5, 415.05, 413.85, 414, 414.75, 415.5, 415.5, 416.1, 417.15, 417.75, 417.45, 417, 416.55, 416.1, 415.8, 415.35, 414.9, 415.5, 415.5, 415.8, 417, 417.45, 417.3, 417.75, 417.45, 417.6, 417.6, 416.7, 416.7, 417.3, 417.6, 417.9, 419.4, 420.3, 420.6, 422.1, 422.85, 422.7, 423.9, 423.75, 423.6, 423.75, 423.3, 423.45, 423.9, 423.15, 423.3, 423.45, 424.2, 424.5, 424.65, 425.1, 426, 426.3, 426.3, 426.15, 425.55, 424.65, 424.95, 424.5, 423.75, 424.2, 423.15, 423.9, 424.35, 423.6, 424.8, 424.95, 423.6, 425.25, 424.8, 423.9, 424.05, 423.3, 423.45, 424.35, 423.75, 422.85, 424.05, 424.35, 423, 422.1, 421.5, 420.6, 420.9, 420.45, 420.3, 421.5, 422.1, 423.3, 424.05, 424.05, 423.9, 424.65, 423.45, 423.45, 422.25, 421.8, 422.55, 422.25, 420.6, 421.05, 420.45, 421.05, 420.75, 419.7, 419.4, 420.15, 419.7, 419.55, 418.5, 418.95, 419.1, 419.1, 418.05, 418.35, 417.75, 417.6, 416.25, 416.1, 416.1, 415.95, 414.9, 414.9, 414.6, 415.05, 414.15, 412.95, 412.65, 413.4, 413.1, 412.95, 412.8, 413.55, 414.15, 414.9, 414.3, 414.15, 414.3, 414.3, 414.15, 413.7, 413.4, 414.15, 414.3, 414.45, 413.4, 412.65, 412.65, 413.1, 412.8, 412.95, 412.5, 413.7, 414.3, 415.65, 415.35, 415.35, 415.2, 416.1, 416.1, 415.95, 416.4, 416.55, 417.3, 417.75, 418.5, 419.1, 419.7, 418.35, 418.8, 418.2, 418.65, 418.05, 416.85, 416.85, 418.5, 418.35, 419.1, 420.3, 420, 420.9, 421.2, 420.6, 420.75, 420.75, 420.45, 421.8, 422.7, 423.3, 423.6, 423.75, 423.15, 423.3, 423.15, 422.85, 422.7, 422.85, 423.45, 424.65, 424.35, 424.2, 425.4, 425.7, 424.95, 424.65, 423.45, 423.3, 423.75, 422.7, 422.85, 423, 422.7, 424.05, 423.9, 423.3, 424.05, 423.9, 423.9, 424.05, 423.75, 423.75, 423.6, 423.15, 423, 423.45, 423.3, 422.1, 422.1, 422.25, 422.1, 421.2, 420.45, 420.6, 421.5, 421.05, 421.05, 420.9, 421.95, 421.35, 420.6, 419.55, 419.7, 419.25, 418.95, 418.65, 418.05, 418.35, 418.35, 418.05, 418.35, 417.9, 417.45, 417.6, 417.75, 417.9, 417.15, 417, 416.25, 415.65, 415.05, 414.9, 414.9, 414.15, 413.55, 414, 413.7, 414.6, 412.8, 413.1, 413.55, 411.75, 411.6, 411.3, 411.6, 412.5, 411.75, 412.05, 413.25, 413.4, 414.3, 414, 414.15, 415.05, 415.2, 415.05]
                xTest = [-165.75, -142.65, -190.2, -281.25, -361.65, -384.75, -341.7, -254.4, -176.85, -149.85, -191.55, -280.95, -365.7, -386.7, -354.3, -263.7, -183.6, -154.95, -200.1, -291.45, -365.55, -398.25, -356.25, -268.5, -187.35, -161.25, -204.15, -298.35, -373.8, -409.5, -367.65, -247.35, -187.5]
                yTest = [256.65, 172.35, 94.05, 64.8, 106.65, 187.35, 265.8, 294.75, 258.15, 173.1, 91.95, 61.2, 117, 200.1, 275.55, 306.75, 267.45, 184.65, 102.75, 75.6, 114, 198.6, 277.8, 305.4, 270.9, 186.45, 103.05, 77.7, 117.75, 197.25, 272.55, 270.75, 266.1]
                zTest = [418.05, 424.8, 421.65, 418.2, 412.8, 410.7, 407.25, 412.65, 421.65, 425.55, 421.35, 410.1, 416.1, 419.1, 419.55, 425.1, 430.5, 435, 433.2, 430.5, 421.65, 418.35, 420.9, 427.35, 432.15, 435.75, 435.45, 429.6, 420.3, 414.45, 422.1, 437.25, 429.15]
                break
*/
            case "yup70": // mounted vertically so Y-axis pointing up; dip=70
                scanTimes = [22061, 22081, 22101, 22121, 22141, 22161, 22181, 22201, 22221, 22241, 22261, 22281, 22301, 22321, 22341, 22361, 22381, 22401, 22423, 22445, 22465, 22485, 22505, 22525, 22545, 22565, 22585, 22605, 22625, 22645, 22665, 22685, 22705, 22725, 22745, 22765, 22785, 22805, 22825, 22845, 22865, 22885, 22905, 22925, 22945, 22965, 22985, 23005, 23025, 23045, 23065, 23085, 23105, 23125, 23145, 23165, 23185, 23205, 23225, 23245, 23265, 23285, 23305, 23325, 23345, 23369, 23389, 23409, 23429, 23449, 23469, 23489, 23509, 23529, 23549, 23569, 23589, 23609, 23629, 23649, 23669, 23689, 23709, 23729, 23749, 23769, 23789, 23809, 23829, 23849, 23869, 23889, 23909, 23929, 23949, 23969, 23989, 24009, 24029, 24052, 24073, 24093, 24113, 24133, 24153, 24173, 24193, 24213, 24233, 24253, 24273, 24293, 24313, 24333, 24353, 24373, 24393, 24413, 24437, 24457, 24477, 24497, 24517, 24537, 24557, 24577, 24597, 24617, 24637, 24657, 24677, 24697, 24717, 24737, 24757, 24777, 24797, 24817, 24837, 24861, 24881, 24901, 24921, 24941, 24961, 24981, 25001, 25021, 25041, 25061, 25081, 25101, 25121, 25141, 25161, 25185, 25205, 25225, 25245, 25265, 25285, 25305, 25325, 25345, 25365, 25385, 25405, 25429, 25449, 25469, 25489, 25509, 25529, 25549, 25569, 25589, 25613, 25637, 25657, 25677, 25697, 25717, 25737, 25757, 25777, 25797, 25817, 25837, 25857, 25877, 25897, 25917, 25937, 25957, 25977, 25997, 26017, 26037, 26062, 26085, 26105, 26125, 26145, 26165, 26185, 26205, 26225, 26245, 26265, 26285, 26305, 26325, 26345, 26365, 26385, 26405, 26425, 26445, 26465, 26485, 26505, 26525, 26545, 26565, 26585, 26605, 26625, 26645, 26665, 26685, 26705, 26725, 26745, 26765, 26785, 26805, 26825, 26845, 26865, 26885, 26905, 26925, 26945, 26965, 26985, 27005, 27025, 27045, 27065, 27085, 27105, 27125, 27145, 27165, 27185, 27212, 27233, 27253, 27273, 27293, 27313, 27333, 27353, 27373, 27393, 27413, 27433, 27453, 27473, 27493, 27513, 27533, 27553, 27573, 27593, 27613, 27633, 27653, 27673, 27693, 27713, 27733, 27753, 27773, 27793, 27813, 27833, 27853, 27873, 27893, 27913, 27933, 27953, 27973, 28001, 28021, 28041, 28061, 28081, 28101, 28121, 28141, 28161, 28181, 28201, 28221, 28241, 28261, 28281, 28301, 28321, 28341, 28361, 28381, 28401, 28421, 28441, 28461, 28481, 28501, 28521, 28548, 28569, 28589, 28609, 28629, 28649, 28669, 28689, 28709, 28729, 28749, 28769, 28789, 28809, 28829, 28849, 28869, 28896, 28917, 28937]
                xData = [-126.3, -124.5, -123.75, -122.7, -121.95, -120.3, -118.2, -115.95, -112.65, -108.45, -104.7, -100.8, -96.45, -92.1, -87.15, -82.35, -78.3, -71.85, -65.25, -59.25, -53.25, -46.95, -41.85, -34.8, -29.7, -24.3, -18.9, -12.9, -7.8, -1.2, 3.6, 9.3, 14.4, 19.5, 23.4, 28.05, 31.65, 36.6, 41.25, 45.3, 49.65, 54.9, 59.7, 64.05, 68.55, 72.3, 76.2, 80.4, 83.85, 86.55, 89.4, 91.5, 93.45, 96, 98.4, 99.45, 101.4, 103.5, 103.8, 105.15, 105.15, 104.55, 106.05, 105.9, 103.8, 103.35, 102.45, 101.7, 101.25, 99.15, 97.35, 96, 93.3, 89.1, 85.2, 80.55, 76.8, 72.6, 68.25, 64.35, 60.45, 56.85, 51.9, 46.8, 42.6, 37.05, 32.85, 27.15, 20.55, 15.45, 8.4, 0, -6.6, -13.8, -19.05, -25.5, -33.3, -39.75, -45.6, -51.45, -58.5, -65.7, -70.95, -75.75, -79.65, -85.05, -90.3, -94.8, -99.75, -104.85, -109.5, -113.25, -115.5, -118.5, -121.65, -123, -124.5, -124.8, -126, -126.75, -126.15, -125.4, -124.35, -123.6, -122.7, -120.9, -119.7, -117.75, -114.9, -112.95, -110.1, -108.3, -105.75, -100.65, -96.15, -92.7, -88.5, -83.7, -78.6, -74.1, -69.15, -64.8, -60, -54.9, -50.25, -44.85, -39.15, -35.55, -30.75, -25.05, -19.35, -13.2, -7.8, -2.1, 3.3, 8.55, 13.8, 18.6, 23.85, 29.4, 33.3, 38.1, 43.35, 47.7, 52.5, 56.85, 61.05, 65.4, 69.9, 73.05, 76.35, 79.05, 82.05, 84.75, 87.9, 90.3, 92.7, 95.1, 97.8, 99.45, 100.65, 101.85, 103.35, 104.85, 105.45, 105.6, 106.2, 107.1, 107.1, 106.2, 105.75, 105.9, 105.75, 103.8, 101.55, 99.9, 98.55, 97.2, 94.35, 91.65, 89.4, 86.7, 84.9, 81.45, 76.2, 73.05, 67.95, 64.2, 59.7, 54.45, 50.4, 47.85, 43.5, 40.5, 35.85, 32.25, 27, 21, 14.4, 8.1, 1.5, -4.65, -10.5, -16.65, -23.1, -28.35, -34.2, -40.95, -47.7, -54.3, -60, -64.5, -69.9, -76.05, -79.5, -83.25, -87.75, -91.65, -96.45, -100.2, -102.3, -105.9, -109.05, -111.75, -114.9, -116.4, -118.95, -121.95, -123.3, -123.45, -125.1, -125.25, -125.85, -125.4, -124.5, -124.05, -124.5, -123, -123, -122.7, -122.55, -121.65, -121.05, -119.85, -117.9, -115.35, -112.2, -109.8, -106.95, -104.1, -101.4, -98.55, -95.4, -93, -88.95, -85.2, -81.3, -77.1, -73.05, -69.75, -65.7, -61.8, -57.9, -53.25, -49.2, -45.75, -40.95, -35.7, -31.2, -26.1, -21.15, -16.2, -10.2, -4.35, 1.5, 7.5, 13.5, 18.6, 24.45, 30, 34.95, 38.55, 41.1, 44.25, 47.85, 51.6, 54.75, 58.05, 61.8, 66.3, 70.2, 73.65, 76.35, 79.5, 82.2, 84.75, 87.15, 88.8, 90.9, 93.6, 94.8, 96.3, 98.25, 100.65, 102.9, 104.55, 104.85, 105.9, 106.65, 107.4, 107.1, 106.05, 105, 104.55, 103.95, 103.2, 101.1, 99.3, 97.2, 94.95, 91.5]
                yData = [-833.85, -833.55, -833.4, -833.55, -833.4, -833.1, -832.8, -832.8, -832.5, -832.5, -831.75, -832.35, -832.35, -832.8, -833.25, -833.1, -832.5, -832.65, -832.2, -831.9, -831.3, -831, -831.3, -831.45, -832.05, -832.05, -832.05, -831.6, -831.6, -831.75, -832.05, -831.75, -831.45, -832.05, -832.2, -832.65, -832.8, -832.65, -833.55, -833.4, -833.25, -833.85, -833.4, -833.7, -833.7, -833.1, -833.55, -833.85, -834.15, -834.45, -834.3, -834.15, -834.15, -835.2, -835.65, -835.35, -835.35, -835.35, -836.1, -836.4, -836.25, -836.4, -836.7, -837, -837.15, -837.6, -837.6, -837.6, -837.45, -838.65, -838.95, -839.85, -840.6, -841.05, -841.35, -842.25, -841.95, -841.95, -841.95, -841.8, -841.8, -841.8, -841.2, -841.05, -841.35, -841.05, -840.6, -841.35, -841.5, -842.4, -842.4, -841.95, -842.4, -842.1, -841.2, -840.45, -840.15, -839.7, -840.15, -839.4, -839.55, -840.3, -840.75, -839.85, -839.7, -839.1, -838.95, -838.2, -837.6, -837, -837, -837.15, -837, -836.55, -836.7, -836.1, -835.95, -835.35, -834.45, -833.1, -833.7, -833.85, -833.55, -832.95, -832.5, -832.8, -833.25, -832.5, -831.3, -831, -831.3, -831.15, -831, -832.35, -832.05, -833.1, -833.85, -834, -833.85, -834, -832.5, -833.1, -832.8, -832.35, -832.2, -832.8, -832.5, -833.4, -832.5, -831.75, -832.05, -832.2, -831.75, -831.45, -830.85, -831.75, -832.05, -832.05, -831.75, -832.35, -833.25, -833.7, -832.8, -833.4, -833.55, -832.8, -832.95, -832.5, -832.35, -832.5, -831.9, -831.75, -832.95, -833.1, -833.1, -833.85, -834.15, -834.6, -834.6, -834.6, -834.3, -834, -833.1, -833.55, -833.7, -834.75, -834.3, -835.5, -837.3, -837.75, -838.35, -839.4, -838.95, -839.4, -839.55, -839.4, -840, -840.6, -840.15, -841.5, -842.7, -842.7, -842.7, -843.3, -842.55, -842.4, -842.1, -841.8, -841.65, -841.8, -842.1, -843, -844.65, -843.3, -844.05, -844.35, -844.65, -844.65, -845.1, -843.3, -843.3, -842.4, -841.5, -840.9, -840, -838.95, -838.65, -839.1, -839.55, -840, -840, -839.85, -839.55, -840, -840.6, -838.95, -837.9, -836.85, -837.15, -837.3, -837.3, -837, -838.05, -837.75, -837.75, -838.35, -837.6, -837.45, -836.4, -835.5, -836.1, -836.1, -834.3, -833.7, -832.95, -833.1, -831.9, -831.6, -831.15, -831, -831, -830.55, -829.65, -830.7, -830.4, -831.45, -830.85, -830.7, -831.6, -831.3, -831, -831, -829.95, -830.85, -831, -830.85, -831.6, -831.75, -832.05, -832.35, -832.95, -832.5, -831.75, -831.3, -830.7, -830.7, -830.4, -829.5, -829.2, -829.05, -829.35, -829.8, -829.95, -829.8, -829.95, -831, -831.6, -831.45, -831, -831.3, -832.2, -831.6, -831.15, -831.15, -832.2, -832.8, -832.2, -831.6, -831.45, -832.05, -832.5, -831.15, -830.85, -831, -831.75, -832.65, -832.95, -832.95, -833.85, -834.75, -835.05, -835.05, -834.9, -835.65, -836.7, -837.45, -837.3, -837.15, -837, -837.6, -837.9, -837.6, -837.6, -837.75, -838.05, -838.2, -839.25, -840.45, -840.15, -840.45, -841.05]
                zData = [-311.7, -306.15, -300.15, -295.5, -289.65, -285.45, -280.65, -274.95, -270, -264.15, -259.2, -255, -249.6, -243.75, -240.6, -236.25, -233.55, -229.95, -226.5, -223.35, -222.3, -219.9, -218.4, -216, -215.55, -215.1, -215.55, -214.8, -214.5, -213.75, -214.2, -214.05, -215.25, -216.3, -217.35, -219.9, -222.45, -225.75, -227.55, -228.75, -230.55, -233.55, -235.2, -239.55, -242.4, -246.6, -250.8, -254.55, -259.35, -263.85, -266.25, -270, -274.65, -279.75, -285.15, -289.35, -294.15, -298.95, -304.65, -309.6, -314.55, -319.2, -323.25, -328.05, -333.3, -338.55, -343.35, -347.55, -352.8, -358.5, -364.65, -370.8, -375.45, -381.45, -387.3, -393.15, -397.95, -401.7, -405.45, -410.7, -414.15, -417.3, -418.8, -421.8, -424.5, -427.05, -427.95, -429, -430.2, -432.3, -433.8, -434.55, -434.4, -434.7, -433.95, -434.1, -433.65, -431.1, -429, -427.5, -424.05, -421.8, -417.9, -414, -409.95, -406.2, -401.55, -396.3, -391.05, -385.65, -379.2, -374.55, -368.55, -362.7, -357.6, -352.2, -347.4, -342.6, -337.2, -330.9, -325.35, -319.8, -314.4, -308.7, -302.4, -296.25, -291.6, -285.45, -279.6, -274.05, -267.6, -263.4, -259.05, -253.8, -250.8, -247.65, -243.75, -240.9, -236.55, -232.8, -229.5, -225.6, -222.45, -220.2, -219, -217.2, -215.55, -214.65, -214.65, -214.65, -214.8, -213.3, -214.05, -213.75, -215.25, -215.7, -216, -216, -217.35, -219.15, -221.7, -222.3, -224.25, -226.95, -230.7, -234.75, -237.9, -241.2, -244.95, -247.95, -250.95, -253.8, -255.75, -258.3, -261.15, -265.35, -269.7, -274.2, -278.85, -284.7, -291.15, -297.75, -302.55, -307.5, -311.85, -317.25, -321.75, -324.45, -327.45, -331.95, -335.55, -340.95, -345, -349.2, -354.3, -358.8, -364.5, -370.05, -373.8, -376.8, -381.15, -385.8, -390.6, -393.9, -397.2, -400.05, -404.4, -408, -411.15, -414.45, -415.8, -417.6, -421.35, -423.6, -424.95, -426.6, -427.65, -430.05, -431.85, -432.15, -432.6, -433.2, -433.8, -434.7, -434.55, -433.95, -433.5, -432.45, -431.85, -429.9, -426.75, -423.3, -420.75, -418.2, -416.4, -412.8, -408.15, -404.7, -401.25, -397.95, -393.3, -387.45, -381.9, -377.25, -373.05, -368.7, -363.6, -358.5, -352.65, -348.3, -343.65, -336.9, -331.35, -326.1, -321.45, -317.25, -312.3, -307.8, -304.35, -300.15, -296.25, -291.75, -287.7, -283.95, -280.05, -275.25, -270.75, -264.6, -260.4, -256.05, -252, -247.8, -244.65, -242.4, -241.65, -238.35, -237.3, -234.75, -232.95, -231.3, -226.95, -224.7, -223.5, -220.65, -218.85, -217.65, -215.7, -216.75, -215.1, -214.65, -214.5, -214.05, -214.5, -215.55, -216, -216.3, -216.75, -217.65, -219.15, -220.35, -220.95, -222.15, -224.85, -227.1, -229.05, -231.75, -233.7, -237, -240.3, -244.05, -247.05, -250.05, -253.8, -258.45, -262.8, -266.55, -269.7, -274.2, -279.45, -282.6, -285.6, -288.3, -292.5, -297.6, -301.8, -305.55, -310.35, -314.7, -320.7, -324.75, -328.2, -331.95, -337.05, -342.75, -348.6, -351.75, -356.55, -361.2, -366.6, -371.55, -375]
                xTest = [98.7, 95.55, 99.45, 27.6, -63.15, -117.9, -115.35, -50.7, 36.9, 95.85, 96.9, 29.25, -63, -120, -118.2, -52.95, 36, 96.45]
                yTest = [-835.8, -832.65, -838.2, -838.8, -838.65, -835.65, -831.15, -828.6, -833.85, -835.2, -840.15, -843.6, -840.15, -832.8, -832.5, -831.9, -833.25, -832.65]
                zTest = [-284.4, -287.85, -370.95, -427.05, -425.55, -360.45, -276.15, -219, -224.25, -282.15, -366.45, -432.75, -421.8, -360.9, -276.3, -216.15, -221.1, -281.55]
                break

 /*           case "tldown70": // mounted at a strange angle with top-left corner downwards; dip=70
                scanTimes = [27893, 27913, 27933, 27953, 27973, 27993, 28013, 28033, 28053, 28073, 28093, 28113, 28133, 28153, 28173, 28193, 28213, 28235, 28257, 28277, 28297, 28317, 28337, 28357, 28377, 28397, 28417, 28437, 28457, 28477, 28497, 28517, 28537, 28557, 28577, 28597, 28617, 28637, 28657, 28677, 28697, 28717, 28737, 28757, 28777, 28797, 28817, 28837, 28857, 28877, 28897, 28917, 28937, 28957, 28977, 28997, 29017, 29037, 29057, 29077, 29097, 29117, 29137, 29157, 29180, 29201, 29221, 29241, 29261, 29281, 29301, 29321, 29341, 29361, 29381, 29401, 29421, 29441, 29461, 29481, 29501, 29521, 29541, 29561, 29581, 29601, 29621, 29641, 29661, 29681, 29701, 29721, 29741, 29761, 29781, 29801, 29821, 29845, 29865, 29885, 29905, 29925, 29945, 29965, 29985, 30005, 30025, 30045, 30065, 30085, 30105, 30125, 30145, 30165, 30185, 30205, 30229, 30249, 30269, 30289, 30309, 30329, 30349, 30369, 30389, 30409, 30429, 30449, 30469, 30489, 30509, 30529, 30549, 30569, 30589, 30609, 30629, 30652, 30673, 30693, 30713, 30733, 30753, 30773, 30793, 30813, 30833, 30853, 30873, 30893, 30913, 30933, 30953, 30973, 30996, 31017, 31037, 31057, 31077, 31097, 31117, 31137, 31157, 31177, 31197, 31217, 31240, 31261, 31281, 31301, 31321, 31341, 31361, 31381, 31401, 31424, 31445, 31469, 31489, 31509, 31529, 31549, 31569, 31589, 31609, 31629, 31649, 31669, 31689, 31709, 31729, 31749, 31769, 31789, 31809, 31829, 31849, 31869, 31897, 31917, 31937, 31957, 31977, 31997, 32017, 32037, 32057, 32077, 32097, 32117, 32137, 32157, 32177, 32197, 32217, 32237, 32257, 32277, 32297, 32317, 32337, 32357, 32377, 32397, 32417, 32437, 32457, 32477, 32497, 32517, 32537, 32557, 32577, 32597, 32617, 32637, 32657, 32677, 32697, 32717, 32737, 32757, 32777, 32797, 32817, 32837, 32857, 32877, 32897, 32917, 32937, 32957, 32977, 32997, 33017, 33044, 33065, 33085, 33105, 33125, 33145, 33165, 33185, 33205, 33225, 33245, 33265, 33285, 33305, 33325, 33345, 33365, 33385, 33405, 33425, 33445, 33465, 33485, 33505, 33525, 33545, 33565, 33585, 33605, 33625, 33645, 33665, 33685, 33705, 33725, 33745, 33765, 33785, 33805, 33833, 33853, 33873, 33893, 33913, 33933, 33953, 33973, 33993, 34013, 34033, 34053, 34073, 34093, 34113, 34133, 34153, 34173, 34193, 34213, 34233, 34253, 34273, 34293, 34313, 34333, 34353, 34380, 34401, 34421, 34441, 34461, 34481, 34501, 34521, 34541, 34561, 34581, 34601, 34621, 34641, 34661, 34681, 34701, 34728, 34749, 34769]
                xData = [-245.55, -243, -239.85, -237.3, -233.85, -230.4, -227.25, -222.45, -217.95, -213.3, -209.1, -204.45, -199.65, -195.9, -192.15, -187.8, -184.65, -178.8, -174.15, -169.35, -163.05, -156.6, -149.4, -142.5, -136.2, -129.75, -122.85, -116.55, -111, -106.2, -100.65, -95.85, -90, -85.35, -80.55, -75.3, -70.65, -66.3, -61.95, -58.05, -53.7, -50.7, -48.45, -45.15, -42.9, -40.35, -38.55, -37.5, -35.55, -34.2, -35.1, -34.8, -34.65, -34.95, -34.95, -35.7, -36.75, -37.35, -39, -41.4, -44.1, -47.7, -51.9, -55.65, -59.7, -64.2, -68.85, -73.65, -77.25, -81, -85.2, -89.85, -93.45, -98.25, -102.9, -108.45, -112.65, -118.05, -123.3, -128.7, -133.65, -139.65, -145.65, -153.3, -159.6, -165, -172.2, -179.85, -186.45, -193.5, -199.5, -205.5, -212.25, -218.1, -222.9, -228, -232.5, -236.85, -241.2, -244.2, -247.5, -250.05, -252.75, -255.45, -257.1, -259.5, -261.45, -262.95, -264.75, -265.2, -265.2, -266.7, -265.65, -265.8, -263.85, -261.45, -258.9, -255.9, -252.45, -249.3, -246, -243, -239.7, -235.65, -232.8, -228.3, -223.8, -218.4, -213.6, -208.5, -204.45, -198.9, -193.35, -188.4, -184.05, -178.2, -173.4, -168.15, -162, -157.2, -151.5, -144.3, -139.05, -132.9, -126.6, -121.05, -115.35, -109.95, -105.6, -100.95, -95.85, -89.7, -84.9, -80.1, -74.55, -70.8, -65.1, -60.75, -58.8, -55.05, -51.45, -49.05, -45.6, -43.8, -41.55, -38.7, -37.2, -35.85, -34.65, -33.6, -33.75, -34.65, -34.8, -35.1, -36, -36.9, -38.25, -39.6, -40.95, -43.8, -46.95, -49.5, -53.25, -57.15, -61.5, -66.45, -70.95, -76.65, -82.05, -88.05, -94.05, -98.1, -102.15, -106.5, -110.55, -114.9, -118.5, -122.4, -127.5, -131.85, -136.2, -141.3, -146.1, -151.2, -156, -161.1, -166.8, -172.5, -177.3, -183.15, -189.15, -193.65, -199.2, -204, -208.05, -212.7, -217.2, -220.65, -225.3, -228.9, -232.5, -237, -240, -242.55, -245.85, -248.7, -251.1, -253.05, -254.7, -256.2, -259.05, -260.4, -261.9, -262.95, -264.15, -265.05, -266.4, -265.8, -266.7, -265.8, -265.5, -264.75, -263.4, -261.9, -260.25, -257.55, -256.5, -253.5, -251.55, -249.15, -245.55, -242.25, -238.5, -234.6, -231.6, -226.95, -223.2, -219.3, -214.8, -211.2, -207, -201.9, -198.15, -193.65, -189.6, -185.55, -180.15, -174, -168.75, -163.05, -157.35, -151.05, -144.15, -138.15, -133.5, -127.95, -122.1, -116.1, -111.6, -107.25, -102.45, -97.05, -92.4, -87.3, -82.05, -75.75, -70.35, -65.1, -60.3, -56.1, -52.2, -49.35, -46.2, -43.35, -41.55, -39.6, -37.35, -35.4, -34.05, -33.75, -33.9, -33, -33.75, -35.7, -36.75, -38.85, -40.35, -41.7, -44.7, -46.5, -49.05, -53.1, -56.1, -59.25, -62.7, -66.3, -70.2, -73.35, -76.95, -80.4, -85.35, -91.05, -97.2, -102.9, -108.9, -113.85, -120.15, -126.15, -133.2, -137.55, -143.7, -148.65, -156, -161.1, -166.05, -169.5, -175.35, -178.8, -183, -187.5, -193.05]
                yData = [353.25, 356.55, 359.55, 362.4, 365.4, 368.85, 371.4, 373.5, 375.45, 376.8, 378.9, 380.85, 381.45, 382.5, 384.15, 385.05, 386.7, 387.45, 387.45, 388.65, 388.8, 387.75, 386.85, 385.65, 383.55, 381.9, 378.9, 376.65, 374.1, 371.85, 369, 366, 363, 360.3, 356.1, 352.35, 347.55, 342.9, 338.4, 332.25, 325.95, 319.5, 312.9, 307.5, 302.25, 296.1, 289.35, 284.1, 279, 274.8, 270, 263.85, 258.3, 253.95, 248.55, 243.15, 236.85, 230.85, 225.45, 220.2, 215.4, 210.3, 205.35, 201.45, 198, 194.25, 190.2, 187.8, 185.4, 182.7, 180, 176.4, 174.75, 173.55, 171.15, 169.35, 168.3, 167.7, 167.25, 165.75, 166.2, 167.4, 168, 169.5, 171.3, 173.55, 175.95, 178.35, 180.3, 183, 186.45, 189.75, 193.8, 199.05, 202.95, 207.6, 212.7, 216.9, 221.1, 226.2, 230.7, 235.5, 240.45, 245.85, 251.7, 258.45, 263.4, 270, 276.75, 283.2, 289.65, 295.8, 301.65, 309.15, 314.7, 320.55, 327.3, 333.6, 339.15, 344.25, 348.75, 353.7, 358.5, 361.8, 364.2, 368.1, 371.85, 373.65, 376.5, 378.6, 380.25, 382.65, 384.45, 384.75, 387.15, 387.6, 387.9, 388.65, 388.8, 387.75, 387.9, 387, 385.8, 384.6, 382.5, 380.85, 378.45, 376.2, 373.05, 370.2, 367.05, 364.35, 360, 356.4, 350.85, 346.2, 341.1, 336.3, 330.3, 325.95, 320.55, 315.75, 310.5, 305.1, 298.95, 293.4, 288, 282.15, 275.4, 269.7, 263.85, 258.15, 252.15, 245.25, 239.1, 235.05, 229.8, 224.4, 219.75, 215.4, 211.5, 207.45, 203.25, 198.6, 194.4, 190.95, 187.35, 183.3, 181.05, 178.05, 175.8, 174.6, 171.9, 169.8, 169.8, 169.05, 168.6, 168.75, 168, 168.9, 170.1, 170.55, 170.1, 171, 172.05, 173.25, 173.85, 174.9, 175.65, 177.6, 178.95, 180.75, 183.45, 186.75, 190.05, 193.95, 198, 201.45, 204.75, 208.05, 211.5, 215.1, 218.1, 221.1, 226.2, 231, 236.1, 241.05, 246.3, 250.95, 256.95, 261.15, 266.25, 271.35, 276.6, 281.25, 286.35, 290.55, 295.2, 299.85, 304.35, 307.95, 313.65, 318.75, 324.45, 329.85, 335.1, 339.9, 345.9, 349.05, 353.4, 357.3, 361.35, 364.2, 367.65, 370.35, 373.5, 376.05, 378.15, 379.65, 381.6, 383.55, 384, 385.8, 386.1, 386.7, 387.75, 388.2, 387.45, 388.35, 387.45, 388.35, 387.9, 386.55, 385.95, 384.45, 382.95, 380.55, 377.1, 374.25, 372, 368.4, 366.15, 361.05, 357, 353.25, 348.75, 343.35, 337.65, 331.05, 326.7, 321.3, 315.3, 309.15, 303.15, 295.35, 288.9, 281.7, 274.35, 268.05, 261.6, 254.85, 250.05, 243.45, 238.05, 232.8, 226.95, 222.15, 217.5, 213, 209.55, 205.2, 201.45, 197.55, 192.9, 188.85, 185.4, 182.25, 180.6, 178.8, 176.25, 175.05, 174.15, 173.55, 172.8, 170.4, 169.2, 169.8, 170.25, 170.25, 170.85, 171.3, 172.5, 174, 174.15, 175.2, 177.15, 179.1, 181.05, 183.75, 186]
                zData = [388.35, 386.25, 383.7, 381.15, 378, 375.75, 374.4, 370.65, 368.4, 365.7, 364.2, 362.1, 358.5, 355.95, 353.55, 350.25, 347.7, 345.45, 343.5, 342.6, 339.45, 338.1, 336.75, 335.25, 332.55, 329.4, 326.7, 326.4, 325.2, 324.45, 323.7, 322.8, 323.4, 324.3, 325.2, 325.35, 324.6, 325.35, 327, 328.65, 329.1, 329.1, 329.55, 330.9, 332.85, 333.45, 334.35, 336.15, 337.8, 340.65, 344.1, 345.75, 348.45, 350.7, 353.1, 356.7, 359.55, 361.35, 365.1, 368.55, 371.55, 374.55, 376.8, 381, 384.6, 387.45, 390.75, 394.65, 398.1, 401.85, 403.2, 406.65, 409.95, 412.8, 415.8, 418.8, 420.75, 424.5, 426.45, 427.35, 429.45, 430.8, 431.85, 433.8, 435.15, 437.25, 438.9, 439.95, 442.05, 444.9, 447.75, 447.6, 448.65, 449.7, 450.75, 450, 448.95, 447.75, 448.8, 447.6, 446.4, 444.9, 444.3, 441.75, 440.25, 438.15, 435.15, 433.8, 432.15, 427.5, 425.7, 421.8, 418.05, 415.8, 412.8, 409.35, 407.85, 404.55, 399.6, 397.05, 393.3, 388.35, 384.15, 380.1, 376.05, 375, 370.8, 367.05, 363.75, 360.3, 357.6, 356.1, 353.25, 350.7, 349.35, 348.15, 346.05, 344.7, 342, 339.45, 337.8, 335.4, 332.55, 331.35, 329.1, 327.45, 326.4, 324.75, 323.85, 324, 322.8, 322.65, 322.95, 322.35, 322.8, 323.7, 324.45, 324.9, 324.75, 324.9, 326.7, 327.9, 328.65, 329.85, 331.95, 334.2, 336.45, 338.85, 340.35, 342.45, 345, 348.15, 351.45, 353.1, 355.95, 360, 363.45, 366, 369.15, 372.45, 377.4, 380.25, 383.4, 386.4, 390, 392.7, 395.25, 398.55, 402.45, 404.7, 408.3, 410.85, 413.25, 414.75, 416.25, 417.6, 420.3, 421.8, 424.35, 426.45, 429.3, 430.65, 433.5, 435.6, 437.4, 437.7, 438.45, 440.55, 442.35, 441.75, 442.2, 443.4, 444.6, 447.3, 448.05, 448.65, 448.8, 448.2, 447.6, 448.35, 447.9, 447, 446.7, 446.4, 446.4, 445.8, 443.7, 441.6, 440.55, 438.6, 438.3, 437.4, 434.4, 432.9, 430.95, 428.1, 426.3, 423, 420.3, 419.7, 417.75, 414.45, 413.7, 410.85, 408.3, 405.15, 401.55, 399.3, 397.35, 393.45, 390, 387.6, 385.5, 382.65, 378, 375.15, 372.3, 369.75, 367.2, 363.6, 361.2, 359.25, 357, 354.15, 351.45, 349.05, 346.65, 343.8, 341.55, 339.3, 336.45, 335.25, 333, 331.65, 330.45, 329.55, 327.6, 327.75, 326.4, 324.9, 324.45, 323.4, 323.1, 322.2, 321.45, 321.9, 321.75, 321.3, 321.9, 323.1, 324.9, 325.35, 327.15, 329.25, 331.5, 333.6, 335.55, 339, 341.85, 343.2, 346.2, 348.45, 351.15, 353.85, 355.95, 358.8, 362.85, 366.15, 370.5, 374.4, 377.25, 378.45, 382.05, 384.3, 386.85, 388.8, 390.75, 393.6, 397.95, 401.85, 405, 408.3, 412.05, 415.35, 418.5, 421.8, 422.85, 425.25, 427.2, 428.4, 430.8, 432.3, 432.9, 435.6, 437.1, 438.45, 440.25, 441, 442.05, 443.4, 443.55]
                xTest = [-156.45, -154.2, -68.85, -24.15, -50.55, -125.25, -210.9, -255.3, -229.35, -148.2, -63.75, -20.7, -47.85, -123.6, -207.3, -249, -224.1, -146.85, -60.75, -22.35, -47.7, -123.6, -208.35, -253.5, -226.8]
                yTest = [388.8, 387.3, 352.5, 275.25, 201.15, 173.25, 204, 284.25, 363.15, 395.55, 352.8, 277.2, 201, 170.25, 202.05, 288.9, 365.7, 394.35, 354.3, 279.6, 205.05, 169.35, 205.05, 285.45, 361.95]
                zTest = [337.8, 338.85, 322.35, 339.9, 382.8, 424.35, 446.25, 426.45, 378.45, 333, 320.85, 337.35, 381.6, 428.7, 449.7, 428.85, 385.5, 339.45, 315.75, 342, 387.75, 429.9, 451.8, 436.65, 386.4]
                break

            case "tldown0": // top-left corner downwards; dip ~= 0 (equatorial)
                scanTimes = [15365, 15385, 15405, 15425, 15445, 15465, 15485, 15505, 15525, 15545, 15565, 15585, 15605, 15625, 15645, 15665, 15685, 15707, 15729, 15749, 15769, 15789, 15809, 15829, 15849, 15869, 15889, 15909, 15929, 15949, 15969, 15989, 16009, 16029, 16049, 16069, 16089, 16109, 16129, 16149, 16169, 16189, 16209, 16229, 16249, 16269, 16289, 16309, 16329, 16349, 16369, 16389, 16409, 16429, 16449, 16469, 16489, 16509, 16529, 16549, 16569, 16589, 16609, 16629, 16649, 16673, 16693, 16713, 16733, 16753, 16773, 16793, 16813, 16833, 16853, 16873, 16893, 16913, 16933, 16953, 16973, 16993, 17013, 17033, 17053, 17073, 17093, 17113, 17133, 17153, 17173, 17193, 17213, 17233, 17253, 17273, 17293, 17313, 17333, 17356, 17377, 17397, 17417, 17437, 17457, 17477, 17497, 17517, 17537, 17557, 17577, 17597, 17617, 17637, 17657, 17677, 17697, 17720, 17741, 17761, 17781, 17801, 17821, 17841, 17861, 17881, 17901, 17921, 17941, 17961, 17981, 18001, 18021, 18041, 18061, 18081, 18101, 18121, 18141, 18165, 18185, 18205, 18225, 18245, 18265, 18285, 18305, 18325, 18345, 18365, 18385, 18405, 18425, 18445, 18465, 18489, 18509, 18529, 18549, 18569, 18589, 18609, 18629, 18649, 18669, 18689, 18709, 18733, 18753, 18773, 18793, 18813, 18833, 18853, 18873, 18893, 18917, 18941, 18961, 18981, 19001, 19021, 19041, 19061, 19081, 19101, 19121, 19141, 19161, 19181, 19201, 19221, 19241, 19261, 19281, 19301, 19321, 19341, 19366, 19389, 19409, 19429, 19449, 19469, 19489, 19509, 19529, 19549, 19569, 19589, 19609, 19629, 19649, 19669, 19689, 19709, 19729, 19749, 19769, 19789, 19809, 19829, 19849, 19869, 19889, 19909, 19929, 19949, 19969, 19989, 20009, 20029, 20049, 20069, 20089, 20109, 20129, 20149, 20169, 20189, 20209, 20229, 20249, 20269, 20289, 20309, 20329, 20349, 20369, 20389, 20409, 20429, 20449, 20469, 20489, 20516, 20537, 20557, 20577, 20597, 20617, 20637, 20657, 20677, 20697, 20717, 20737, 20757, 20777, 20797, 20817, 20837, 20857, 20877, 20897, 20917, 20937, 20957, 20977, 20997, 21017, 21037, 21057, 21077, 21097, 21117, 21137, 21157, 21177, 21197, 21217, 21237, 21257, 21277, 21305, 21325, 21345, 21365, 21385, 21405, 21425, 21445, 21465, 21485, 21505, 21525, 21545, 21565, 21585, 21605, 21625, 21645, 21665, 21685, 21705, 21725, 21745, 21765, 21785, 21805, 21825, 21852, 21873, 21893, 21913, 21933, 21953, 21973, 21993, 22013, 22033, 22053, 22073, 22093, 22113, 22133, 22153, 22173, 22200, 22221, 22241]
                xData = [20.7, 14.7, 7.2, 0.15, -8.55, -17.25, -26.55, -36.6, -48, -58.5, -69.3, -82.2, -94.8, -108.3, -120.75, -133.65, -147.75, -161.85, -176.25, -192.9, -208.35, -226.65, -244.35, -262.35, -280.35, -296.4, -311.4, -327.15, -341.25, -355.05, -368.55, -381.9, -395.7, -409.05, -422.7, -436.05, -447.6, -458.7, -468.6, -477.6, -485.25, -492, -496.95, -502.8, -506.4, -509.85, -511.95, -513.3, -512.7, -511.95, -509.4, -507.6, -504.6, -501.75, -497.4, -491.7, -484.35, -477.15, -469.35, -460.2, -448.95, -438, -428.55, -419.4, -409.8, -398.1, -386.85, -376.95, -366.45, -354.75, -344.7, -332.85, -322.2, -311.4, -299.7, -288.15, -277.5, -265.05, -252.6, -240.45, -228.45, -216.15, -204.3, -191.4, -177.45, -163.05, -149.1, -135.3, -121.35, -108.15, -95.4, -84.75, -75.45, -65.4, -56.25, -47.55, -38.1, -28.8, -19.65, -10.95, -3, 4.65, 12.45, 18.9, 25.05, 30.75, 35.25, 38.1, 42, 45.15, 47.1, 48.3, 48.75, 47.7, 47.25, 45, 41.85, 38.1, 33, 28.8, 24.45, 18.6, 11.7, 3.6, -5.1, -14.4, -25.95, -37.5, -48.15, -60.6, -73.8, -87, -101.25, -114.9, -130.5, -146.85, -161.85, -178.5, -194.1, -210.15, -229.05, -247.2, -265.65, -284.55, -301.5, -320.1, -337.35, -351.6, -366, -379.8, -393.15, -406.95, -419.7, -432, -444.9, -455.7, -465.9, -475.05, -483.45, -489.6, -496.5, -501.45, -506.25, -510, -513.3, -514.35, -516.15, -515.25, -514.35, -512.7, -510.3, -506.55, -503.1, -498.3, -492.6, -486.3, -479.25, -470.1, -460.5, -449.55, -439.05, -428.55, -417.15, -405.9, -396.15, -386.85, -377.25, -366.6, -356.55, -346.5, -335.25, -323.85, -312.15, -301.95, -291.6, -279.75, -267.9, -256.2, -244.05, -230.55, -215.85, -200.4, -186.6, -173.55, -160.05, -146.25, -133.65, -121.5, -112.05, -101.85, -90.6, -82.05, -74.1, -66.3, -59.1, -50.7, -43.65, -36.9, -29.85, -23.55, -17.55, -11.85, -6.3, 0, 5.1, 10.95, 16.65, 21.45, 26.1, 30.45, 33.75, 37.2, 40.65, 42.6, 45.45, 47.1, 47.55, 47.55, 47.7, 45.75, 43.65, 41.7, 38.85, 35.55, 31.5, 25.8, 20.25, 14.1, 4.95, -3.75, -12.6, -23.55, -33.3, -44.85, -56.85, -69.3, -84.45, -101.4, -116.85, -135.15, -151.8, -169.05, -186.3, -200.85, -215.1, -228.45, -241.8, -254.7, -267.9, -279.75, -292.8, -304.65, -317.25, -329.55, -343.2, -355.95, -369, -382.35, -394.8, -407.25, -419.25, -429.6, -440.4, -450.45, -458.1, -466.65, -475.35, -482.7, -490.05, -495.45, -500.55, -504.9, -509.1, -510.75, -513, -514.2, -514.5, -515.25, -516.15, -515.1, -514.2, -511.8, -510, -507.9, -504.15, -499.35, -494.25, -488.85, -483.6, -476.55, -469.95, -462.9, -456.6, -449.4, -441, -431.25, -421.65, -411.3, -401.7, -390.75, -380.55, -370.65, -361.05, -352.35, -342.6, -331.95, -322.8, -313.8, -304.35, -295.35, -285, -274.8, -265.8, -254.55, -242.7, -229.65, -214.8, -200.55, -186.45, -171.75, -157.5]
                yData = [36.9, 27, 17.4, 7.2, -0.75, -9.6, -17.85, -25.05, -31.95, -37.05, -42.75, -47.85, -52.35, -56.1, -60.15, -61.8, -64.35, -65.1, -66.45, -65.7, -63.3, -60.9, -58.95, -55.35, -50.85, -45.3, -38.85, -33.15, -25.95, -16.95, -7.8, 2.1, 12.75, 24.75, 37.5, 51.45, 64.8, 79.5, 93.75, 108.9, 123.45, 139.2, 154.35, 170.1, 185.4, 200.7, 216.15, 230.55, 244.35, 258.45, 271.95, 284.55, 298.65, 311.25, 322.8, 335.55, 348, 359.85, 371.85, 382.65, 392.7, 403.65, 412.35, 419.55, 427.2, 432.9, 438.45, 443.25, 447.6, 451.35, 455.1, 457.2, 460.65, 462.15, 464.4, 464.85, 464.4, 464.25, 463.2, 461.4, 459.9, 457.05, 453.9, 451.05, 446.4, 441.45, 435.45, 428.55, 421.5, 415.65, 406.95, 400.35, 391.8, 383.25, 374.7, 364.8, 353.7, 343.8, 331.95, 321.45, 310.2, 298.95, 287.55, 275.4, 263.4, 250.65, 238.05, 224.7, 210.45, 197.55, 185.1, 172.35, 160.5, 146.25, 132.3, 118.65, 105.15, 91.2, 77.4, 63.3, 51.45, 39.75, 28.95, 18.15, 7.5, -1.65, -10.95, -19.35, -27.15, -34.5, -41.25, -46.5, -52.35, -55.8, -59.25, -61.8, -63.75, -65.1, -66.3, -65.25, -64.05, -60.75, -55.95, -50.85, -45, -37.05, -29.55, -21.6, -13.35, -4.65, 6.15, 16.95, 27.75, 40.35, 53.7, 67.05, 81.45, 94.8, 110.25, 124.8, 139.35, 153, 168, 183.3, 199.2, 213.6, 229.2, 243.15, 259.2, 274.05, 287.7, 301.2, 314.1, 327.6, 341.55, 353.55, 364.05, 375.45, 385.2, 395.4, 403.95, 412.05, 419.25, 426.45, 431.85, 437.55, 442.8, 446.7, 449.55, 451.95, 454.65, 456, 456.9, 456.75, 456.9, 457.95, 458.55, 457.5, 456.75, 454.95, 451.8, 448.5, 443.25, 438.45, 433.65, 427.95, 422.4, 417.9, 412.35, 406.8, 400.8, 394.05, 387.6, 380.1, 373.2, 366.45, 359.1, 351.3, 343.65, 335.85, 328.65, 319.65, 310.5, 301.5, 293.25, 284.4, 275.25, 264.9, 255.45, 245.7, 235.95, 224.55, 214.05, 203.4, 191.55, 178.8, 166.8, 153.9, 140.55, 126.45, 112.95, 101.1, 89.4, 76.2, 64.5, 53.1, 41.85, 30.6, 19.5, 8.7, -1.5, -11.1, -20.1, -28.5, -36.75, -44.4, -51.45, -56.85, -61.95, -64.95, -67.2, -68.7, -68.85, -68.85, -68.25, -66, -64.95, -62.4, -59.85, -56.7, -51.75, -46.8, -42.3, -36.3, -29.25, -21.6, -14.85, -6.75, 2.55, 11.4, 21.75, 31.35, 42, 54.45, 66.3, 77.55, 90.15, 102.15, 114.75, 127.8, 140.25, 152.7, 166.8, 181.2, 195.9, 210.15, 223.8, 238.2, 251.1, 263.85, 275.7, 286.8, 297.9, 308.85, 318.15, 329.55, 339.9, 348.45, 357.75, 366.45, 374.85, 383.85, 390.3, 396.6, 404.1, 410.7, 417.45, 422.7, 427.65, 433.65, 438.3, 441.6, 444.45, 446.55, 449.55, 452.25, 453.9, 454.95, 456.15, 456.45, 457.8, 458.25, 457.95, 456.9, 455.1, 453.6, 451.05, 447.15, 442.35, 437.1, 431.25]
                zData = [172.35, 181.35, 189.3, 197.55, 206.55, 215.1, 224.1, 232.2, 240.15, 248.4, 255.6, 262.8, 271.2, 278.25, 286.05, 293.25, 300.3, 307.05, 314.1, 320.55, 327.15, 333.6, 340.05, 344.85, 348.9, 352.65, 355.5, 357.45, 359.55, 360.6, 362.1, 363.6, 363.75, 363, 362.4, 360.15, 358.5, 355.5, 352.05, 347.7, 344.25, 338.85, 333.15, 327, 320.55, 314.25, 308.1, 300.3, 292.2, 284.55, 276.75, 269.4, 262.2, 253.35, 245.7, 237.75, 228.45, 218.25, 208.65, 197.85, 189, 179.1, 169.05, 161.25, 154.2, 146.1, 138.6, 131.4, 124.5, 119.7, 114, 107.4, 102.3, 97.2, 91.95, 86.4, 81, 76.65, 72.9, 67.8, 64.2, 60.45, 57.3, 53.55, 50.1, 46.05, 42.75, 39.15, 36.9, 35.25, 33.75, 32.25, 32.85, 33.75, 35.4, 36.75, 37.95, 39.15, 41.1, 43.05, 45.15, 47.25, 48.6, 52.05, 57.3, 60.75, 64.2, 68.7, 73.5, 80.4, 85.5, 90.3, 97.2, 105.6, 114.9, 123.9, 132.3, 140.7, 148.95, 157.8, 165.75, 174, 182.7, 190.65, 199.65, 208.35, 217.95, 227.55, 236.55, 244.8, 254.7, 264.15, 274.5, 283.2, 292.5, 300.6, 309.75, 316.2, 323.25, 330.3, 336, 340.5, 345.3, 349.05, 354.15, 358.05, 359.1, 361.05, 363.3, 364.05, 364.8, 364.2, 363.45, 364.5, 363.15, 360.9, 358.2, 353.55, 349.5, 344.4, 338.1, 333.45, 326.7, 320.7, 315.15, 308.4, 301.65, 294.15, 286.05, 278.25, 269.7, 261.6, 252.45, 243.15, 234.9, 225.15, 215.7, 205.5, 194.85, 185.55, 176.4, 167.7, 159.75, 152.1, 145.65, 139.35, 132.9, 127.05, 120.75, 115.05, 109.5, 103.95, 98.4, 93.15, 88.35, 82.95, 78.3, 74.7, 70.65, 66.3, 63.6, 60, 56.85, 52.95, 48.45, 45.45, 43.8, 41.25, 39.75, 38.4, 37.5, 37.2, 37.65, 36.9, 36.9, 37.2, 37.65, 38.7, 39.9, 41.1, 42.9, 43.8, 45.9, 48, 50.55, 52.65, 54.45, 58.05, 61.65, 64.2, 67.95, 71.7, 75.3, 79.8, 84.6, 91.2, 97.5, 103.05, 109.35, 116.25, 124.65, 131.55, 137.4, 144.45, 154.05, 162.75, 172.35, 181.2, 189.75, 198.45, 208.35, 217.95, 226.2, 235.2, 243.6, 253.2, 264.15, 274.05, 282.3, 292.2, 300.3, 307.35, 315.45, 320.7, 326.4, 331.2, 335.1, 340.2, 344.4, 347.7, 351.9, 353.7, 356.7, 359.55, 360, 362.4, 363.3, 363.6, 364.35, 364.5, 364.35, 364.65, 362.85, 360.75, 357.6, 355.5, 352.8, 350.4, 346.65, 343.35, 339, 336.15, 330.45, 324.45, 316.8, 310.35, 303.45, 297.15, 290.1, 283.5, 276.45, 268.95, 262.05, 254.4, 247.05, 239.25, 231.6, 223.8, 218.25, 212.55, 205.2, 199.2, 192.45, 186, 178.5, 170.55, 162, 155.7, 148.65, 142.35, 135.15, 130.35, 125.4, 120.15, 114.45, 108.6, 104.1, 99.3, 94.5, 89.7, 84.6, 80.7, 76.95, 71.55, 68.4, 63.3, 58.35, 54.75, 51.15, 48.75, 46.65]
                xTest = [-117, -116.85, -342.6, -494.85, -492.9, -328.5, -119.1, 26.7, 31.35, -113.25, -342.15, -502.5, -492, -337.8, -119.4, 25.65, 36, -112.5]
                yTest = [-62.85, -58.8, -27.9, 141.6, 343.05, 458.4, 411.6, 255.45, 56.7, -53.1, -24.15, 148.2, 347.55, 456.9, 419.25, 262.2, 64.2, -54]
                zTest = [282.45, 287.25, 363.75, 340.8, 230.25, 109.5, 41.1, 53.1, 156.15, 281.25, 359.4, 333.15, 224.85, 104.85, 34.35, 51.45, 150.3, 276.9]
                break

            case "zdown15": // mounted horizontally; dip ~= 15
                scanTimes = [4793, 4813, 4833, 4853, 4873, 4893, 4913, 4933, 4953, 4973, 4993, 5013, 5033, 5053, 5073, 5093, 5113, 5137, 5157, 5177, 5197, 5217, 5237, 5257, 5277, 5297, 5317, 5337, 5357, 5377, 5397, 5417, 5437, 5457, 5477, 5497, 5517, 5537, 5557, 5577, 5597, 5617, 5637, 5657, 5677, 5697, 5717, 5737, 5757, 5777, 5797, 5817, 5837, 5857, 5877, 5897, 5917, 5937, 5957, 5977, 5997, 6017, 6037, 6057, 6080, 6101, 6121, 6141, 6161, 6181, 6201, 6221, 6241, 6261, 6281, 6301, 6321, 6341, 6361, 6381, 6401, 6421, 6441, 6461, 6481, 6501, 6521, 6541, 6561, 6581, 6601, 6621, 6641, 6661, 6681, 6701, 6724, 6745, 6765, 6785, 6805, 6825, 6845, 6865, 6885, 6905, 6925, 6945, 6965, 6985, 7005, 7025, 7045, 7065, 7085, 7105, 7128, 7149, 7169, 7189, 7209, 7229, 7249, 7269, 7289, 7309, 7329, 7349, 7369, 7389, 7409, 7429, 7449, 7469, 7489, 7509, 7529, 7549, 7569, 7593, 7613, 7633, 7653, 7673, 7693, 7713, 7733, 7753, 7773, 7793, 7813, 7833, 7853, 7873, 7893, 7916, 7937, 7957, 7977, 7997, 8017, 8037, 8057, 8077, 8097, 8117, 8137, 8157, 8181, 8201, 8221, 8241, 8261, 8281, 8301, 8321, 8341, 8365, 8385, 8405, 8425, 8445, 8465, 8485, 8505, 8525, 8545, 8565, 8585, 8605, 8625, 8645, 8665, 8685, 8705, 8725, 8745, 8765, 8790, 8813, 8833, 8853, 8873, 8893, 8913, 8933, 8953, 8973, 8993, 9013, 9033, 9053, 9073, 9093, 9113, 9133, 9153, 9173, 9193, 9213, 9233, 9253, 9273, 9293, 9313, 9333, 9353, 9373, 9393, 9413, 9433, 9453, 9473, 9493, 9513, 9533, 9553, 9573, 9593, 9613, 9633, 9653, 9673, 9693, 9713, 9733, 9753, 9773, 9793, 9813, 9833, 9853, 9873, 9893, 9913, 9933, 9961, 9981, 10001, 10021, 10041, 10061, 10081, 10101, 10121, 10141, 10161, 10181, 10201, 10221, 10241, 10261, 10281, 10301, 10321, 10341, 10361, 10381, 10401, 10421, 10441, 10461, 10481, 10501, 10521, 10541, 10561, 10581, 10601, 10621, 10641, 10661, 10681, 10701, 10729, 10749, 10769, 10789, 10809, 10829, 10849, 10869, 10889, 10909, 10929, 10949, 10969, 10989, 11009, 11029, 11049, 11069, 11089, 11109, 11129, 11149, 11169, 11189, 11209, 11229, 11249, 11276, 11297, 11317, 11337, 11357, 11377, 11397, 11417, 11437, 11457, 11477, 11497, 11517, 11537, 11557, 11577, 11597, 11624, 11645, 11665]
                xData = [69.45, 77.1, 83.4, 88.5, 92.25, 95.4, 96.15, 95.55, 94.8, 93.45, 90.6, 85.8, 79.65, 72.6, 63.75, 53.85, 43.05, 31.95, 20.4, 7.8, -9.15, -24.3, -41.25, -59.4, -78.45, -98.1, -116.85, -132.9, -150.3, -167.1, -183.75, -200.1, -216, -233.1, -249.45, -266.25, -283.5, -300.75, -318.9, -336.75, -353.25, -369, -382.65, -395.55, -407.4, -417.15, -426, -434.85, -443.1, -452.1, -460.05, -467.55, -474.75, -481.65, -487.05, -493.5, -498, -502.35, -506.7, -510.3, -512.4, -516.3, -518.1, -520.35, -521.55, -522.15, -522.6, -522.75, -521.4, -520.95, -521.25, -519.9, -518.4, -516.45, -515.25, -513.45, -509.7, -504.9, -501.15, -497.1, -492.9, -487.35, -481.95, -476.4, -471.3, -465.6, -459, -451.8, -442.95, -433.65, -423.6, -412.8, -401.55, -390, -378.45, -367.8, -357, -346.35, -333.75, -321.45, -308.7, -295.2, -281.55, -266.4, -249.9, -236.1, -221.1, -207, -192.9, -179.1, -165.9, -155.1, -142.35, -129.75, -116.1, -101.7, -88.95, -74.7, -60.15, -45.9, -32.85, -20.25, -9, 2.85, 12.3, 22.8, 31.5, 41.1, 48.75, 57, 64.5, 70.5, 75, 79.8, 83.85, 87.75, 89.85, 91.5, 93.45, 94.2, 94.35, 91.2, 88.95, 85.95, 81.9, 77.7, 71.7, 66, 60.6, 53.85, 45.9, 37.65, 28.65, 20.1, 10.2, 0.6, -10.5, -21.45, -32.7, -45, -57.45, -70.05, -82.8, -94.5, -105.45, -117.3, -129.45, -140.55, -154.5, -168.6, -182.55, -197.4, -212.25, -225.3, -240.15, -251.55, -263.4, -275.4, -286.8, -297.45, -307.65, -317.4, -328.05, -338.1, -346.95, -356.1, -365.85, -376.5, -385.95, -396, -405.3, -415.65, -425.85, -435.6, -445.8, -455.55, -463.8, -471.9, -479.55, -486, -491.4, -496.2, -500.25, -504.9, -508.65, -512.4, -516.6, -520.05, -522.15, -523.65, -524.4, -524.55, -523.8, -521.85, -519.3, -516.75, -514.5, -511.05, -507.9, -504.6, -500.25, -495.6, -491.25, -485.55, -480.15, -474.6, -466.8, -459.75, -452.1, -443.1, -434.7, -424.05, -411.9, -400.8, -389.1, -376.35, -363.45, -349.8, -337.2, -324.15, -310.65, -297.45, -282.9, -269.1, -254.55, -239.1, -224.85, -210.3, -195, -180.45, -165.3, -151.05, -135, -117.75, -101.55, -85.5, -70.65, -55.2, -39.6, -26.85, -16.05, -4.95, 6.45, 16.95, 27, 36.15, 45.15, 54.6, 63.15, 70.2, 76.05, 82.5, 87.3, 91.8, 95.4, 97.2, 98.4, 99.6, 98.7, 96.9, 93.75, 88.95, 85.95, 82.2, 76.95, 71.25, 64.95, 58.65, 52.5, 43.95, 34.35, 25.35, 14.85, 3.3, -7.8, -19.2, -34.2, -49.05, -64.2, -79.5, -93.3, -109.95, -126.3, -140.1, -154.95, -170.25, -185.4, -201, -215.7, -231.45, -246.15, -260.25, -274.65, -289.5, -304.05, -318.45, -331.2, -344.25, -357, -368.55, -378.45, -388.35, -399.45, -410.4, -421.2, -431.1, -440.85, -450.45, -458.7, -464.85, -471.6, -476.7, -481.35, -485.4, -489.3, -494.4, -498.9, -501.45, -502.95, -506.25, -509.25]
                yData = [335.4, 318, 301.05, 283.5, 266.85, 250.2, 234.3, 218.25, 201.6, 184.8, 167.55, 151.05, 134.85, 117.9, 101.1, 84.3, 68.1, 53.1, 39.45, 25.05, 10.8, -2.4, -14.1, -24.75, -33.75, -43.8, -52.2, -57.9, -63.45, -67.95, -72.15, -75, -75.6, -75.6, -74.55, -71.25, -69.15, -65.1, -60.15, -55.05, -49.65, -43.2, -37.5, -29.1, -21.3, -14.1, -7.35, 0.3, 7.95, 16.05, 23.7, 31.35, 39.9, 49.35, 59.4, 68.7, 78.15, 88.2, 98.85, 109.35, 120.15, 130.65, 141.45, 151.95, 160.8, 170.25, 178.95, 186.9, 194.25, 202.5, 211.65, 221.55, 230.7, 240.15, 250.2, 260.55, 270, 279.3, 288.75, 297.45, 306.75, 315.6, 324, 332.4, 341.4, 350.25, 359.7, 368.1, 378.3, 389.1, 398.85, 407.1, 416.4, 424.35, 431.1, 437.25, 442.2, 447.9, 453.6, 456.75, 461.4, 466.05, 469.5, 473.25, 475.8, 477.6, 480, 480.9, 481.95, 481.5, 480.3, 478.95, 476.85, 474.6, 471.9, 468.45, 465.6, 460.65, 454.65, 448.05, 440.25, 432.3, 423.45, 413.1, 403.65, 393.75, 384.3, 373.95, 361.2, 349.95, 338.7, 325.35, 313.2, 300.3, 287.25, 274.8, 260.85, 246.45, 233.1, 218.55, 201.6, 186, 171.3, 156.3, 140.85, 125.7, 110.25, 97.5, 83.85, 71.1, 58.95, 47.1, 36.15, 26.1, 15.75, 5.7, -5.1, -15.45, -25.95, -35.25, -44.1, -52.65, -60.3, -66.45, -72, -75.75, -80.55, -85.05, -88.5, -92.1, -95.25, -98.1, -99.9, -100.5, -100.95, -100.95, -99.15, -97.35, -95.25, -92.85, -90.75, -87.6, -83.7, -80.55, -77.85, -74.55, -70.65, -65.85, -60.45, -55.35, -49.35, -42.75, -36.45, -28.35, -18.9, -11.4, -2.55, 6.15, 16.5, 27.15, 35.4, 43.2, 52.5, 62.1, 71.85, 82.2, 93.3, 104.4, 115.8, 128.7, 140.85, 153.75, 165.75, 177.9, 190.05, 203.1, 214.95, 226.95, 238.65, 250.05, 261.15, 270.75, 280.65, 290.55, 300.15, 308.85, 317.25, 326.7, 337.35, 346.2, 355.2, 364.2, 374.1, 382.65, 389.7, 396.9, 403.95, 409.8, 415.5, 420.3, 424.8, 429.75, 433.95, 437.4, 441.45, 444.3, 446.4, 447.75, 448.8, 448.65, 448.35, 447.15, 444.3, 440.85, 437.7, 432.15, 426.3, 419.25, 411.6, 404.55, 396.3, 387.3, 379.05, 370.05, 360.45, 349.65, 337.5, 326.1, 313.65, 300.45, 286.95, 273.6, 259.5, 246, 229.5, 214.05, 197.7, 182.1, 166.05, 150.3, 134.4, 120.6, 107.7, 95.85, 83.25, 70.95, 59.4, 48, 37.65, 25.35, 13.95, 3.15, -7.95, -18, -26.85, -36.75, -44.7, -53.55, -61.2, -66.6, -73.65, -78.45, -82.65, -87, -88.65, -90.3, -92.85, -92.1, -92.85, -91.5, -89.7, -87.9, -85.35, -81.3, -77.1, -72.75, -69.45, -63.75, -59.25, -54, -48.15, -42.75, -34.5, -24.45, -15.6, -7.05, 1.5, 9.6, 18.3, 24.9, 31.5, 38.4, 46.35, 54.75, 63.6, 72.6, 81.9, 91.05, 100.8, 110.7, 122.4]
                zData = [99, 99.9, 100.35, 101.4, 102.15, 102, 102, 101.7, 100.95, 100.5, 98.85, 97.2, 96.45, 95.4, 94.65, 93.6, 92.85, 92.85, 93, 91.5, 89.85, 87.75, 84.45, 81.15, 78, 74.25, 72, 69.45, 66.3, 64.95, 62.85, 60, 59.1, 57.6, 55.65, 54.9, 53.4, 52.05, 51, 48.75, 46.65, 45, 42.9, 41.1, 40.2, 39.75, 39.3, 38.7, 39.3, 39.15, 38.55, 37.5, 37.35, 36, 34.8, 33.3, 32.7, 32.25, 32.1, 30.6, 30.3, 29.7, 29.25, 28.35, 27.6, 27.75, 26.7, 26.55, 27.45, 27, 26.25, 26.55, 26.1, 27.6, 27.15, 26.85, 27.6, 28.5, 29.4, 30.3, 31.05, 31.8, 32.55, 32.55, 34.05, 33.9, 33.9, 33.9, 34.8, 35.55, 36.75, 36, 36.9, 36.6, 36.45, 36.6, 36.6, 36.9, 38.85, 39.6, 41.7, 42.75, 44.1, 46.05, 47.25, 47.4, 48, 48, 49.5, 50.4, 50.1, 51.45, 52.35, 53.7, 56.1, 56.85, 58.2, 58.95, 60, 61.95, 63.15, 63.6, 64.35, 65.55, 67.05, 67.5, 68.1, 68.1, 68.55, 70.05, 70.35, 71.55, 71.25, 71.25, 72.15, 71.7, 70.8, 70.2, 69.75, 69.9, 68.7, 67.95, 67.5, 67.05, 66.3, 64.5, 63.9, 63.75, 63, 62.25, 61.65, 61.8, 62.55, 61.95, 61.2, 61.8, 62.7, 62.7, 61.8, 60.75, 60.45, 61.05, 59.7, 58.65, 58.05, 57.3, 56.85, 57.3, 56.4, 56.25, 56.4, 55.8, 54.75, 53.7, 51.75, 49.95, 48.45, 46.65, 46.05, 45.45, 45.15, 44.1, 44.25, 43.8, 43.05, 42, 42.15, 41.1, 40.95, 39.6, 39.45, 39, 38.55, 37.95, 37.35, 36.75, 36.45, 36.15, 36, 35.55, 34.95, 35.7, 34.8, 33.9, 32.4, 31.8, 30.6, 30, 29.1, 29.4, 29.85, 30.45, 30.6, 32.25, 33.15, 33.45, 34.65, 35.55, 36.75, 37.05, 37.35, 38.55, 40.2, 40.65, 42, 43.2, 44.7, 45.6, 46.05, 46.2, 46.65, 48.15, 48.45, 49.65, 50.1, 51.3, 52.95, 55.35, 54.75, 55.95, 57, 58.65, 60.15, 60.45, 60.75, 62.85, 64.05, 64.8, 66, 66.6, 68.7, 69.75, 70.35, 70.65, 72.9, 74.85, 76.2, 78.15, 80.1, 82.8, 85.2, 85.95, 87.75, 89.7, 90.45, 90.9, 91.5, 92.4, 93.6, 93.6, 93.9, 94.95, 95.85, 96.45, 97.5, 98.7, 99.45, 101.7, 102.6, 104.55, 105.9, 108.15, 109.2, 111, 111.45, 112.05, 112.8, 113.4, 113.4, 114.3, 114.75, 115.2, 116.25, 116.4, 117.15, 117.3, 117.6, 117.6, 118.2, 118.2, 118.65, 118.65, 117.75, 115.95, 115.65, 114, 112.05, 110.25, 107.85, 106.35, 105.45, 102.6, 101.1, 100.35, 99.3, 98.85, 98.55, 96.9, 97.05, 96.45, 95.1, 93.9, 92.55, 90.9, 90.6, 89.55, 87.75, 87.3, 87, 87.45, 86.7, 86.4, 85.65, 86.25, 85.95, 85.95, 85.35, 85.8, 85.35, 85.2, 83.85]
                xTest = [41.25, 35.1, -152.25, -372.75, -509.1, -497.25, -307.5, -66.75, 85.5, 38.85, -152.7, -365.7, -498.6, -492.3, -304.35, -66.6, 78.45, 31.8]
                yTest = [62.85, 61.8, -66.45, -37.5, 104.85, 294, 436.35, 418.8, 256.2, 22.95, -103.8, -67.35, 84.6, 301.65, 458.85, 460.2, 290.85, 42.15]
                zTest = [81, 91.8, 67.5, 43.8, 34.8, 28.95, 37.35, 70.35, 81.15, 96.3, 89.55, 82.35, 76.2, 78, 73.35, 73.95, 55.35, 66.45]
                break


*/ 
            case "dashboard70": // mounted as 45-degree tilted dashboard; dip ~= 70 
                // In this instance an adjacent speaker magnet was dominating the Earth's field by
                // a factor af ~150x, yet it is still possible to extract the small Spin-Circle field-vector!
                scanTimes = [16533, 16553, 16573, 16593, 16613, 16633, 16653, 16673, 16693, 16713, 16733, 16753, 16773, 16793, 16813, 16833, 16853, 16873, 16897, 16917, 16937, 16957, 16977, 16997, 17017, 17037, 17057, 17077, 17097, 17117, 17137, 17157, 17177, 17197, 17217, 17237, 17257, 17277, 17297, 17317, 17337, 17357, 17377, 17397, 17417, 17437, 17457, 17477, 17497, 17517, 17537, 17557, 17577, 17597, 17617, 17637, 17657, 17677, 17697, 17717, 17737, 17757, 17777, 17797, 17817, 17841, 17861, 17881, 17901, 17921, 17941, 17961, 17981, 18001, 18021, 18041, 18061, 18081, 18101, 18121, 18141, 18161, 18181, 18201, 18221, 18241, 18261, 18281, 18301, 18321, 18341, 18361, 18381, 18401, 18421, 18441, 18461, 18481, 18501, 18524, 18545, 18565, 18585, 18605, 18625, 18645, 18665, 18685, 18705, 18725, 18745, 18765, 18785, 18805, 18825, 18845, 18865, 18885, 18908, 18929, 18949, 18969, 18989, 19009, 19029, 19049, 19069, 19089, 19109, 19129, 19149, 19169, 19189, 19209, 19229, 19249, 19269, 19289, 19309, 19329, 19352, 19373, 19393, 19413, 19433, 19453, 19473, 19493, 19513, 19533, 19553, 19573, 19593, 19613, 19633, 19653, 19676, 19697, 19717, 19737, 19757, 19777, 19797, 19817, 19837, 19857, 19877, 19897, 19920, 19941, 19961, 19981, 20001, 20021, 20041, 20061, 20081, 20105, 20125, 20145, 20165, 20185, 20205, 20225, 20245, 20265, 20285, 20305, 20325, 20345, 20365, 20385, 20405, 20425, 20445, 20465, 20485, 20505, 20530, 20553, 20573, 20593, 20613, 20633, 20653, 20673, 20693, 20713, 20733, 20753, 20773, 20793, 20813, 20833, 20853, 20873, 20893, 20913, 20933, 20953, 20973, 20993, 21013, 21033, 21053, 21073, 21093, 21113, 21133, 21153, 21173, 21193, 21213, 21233, 21253, 21273, 21293, 21313, 21333, 21353, 21373, 21393, 21413, 21433, 21453, 21473, 21493, 21513, 21533, 21553, 21573, 21593, 21613, 21633, 21653, 21673, 21700, 21721, 21741, 21761, 21781, 21801, 21821, 21841, 21861, 21881, 21901, 21921, 21941, 21961, 21981, 22001, 22021, 22041, 22061, 22081, 22101, 22121, 22141, 22161, 22181, 22201, 22221, 22241, 22261, 22281, 22301, 22321, 22341, 22361, 22381, 22401, 22421, 22441, 22469, 22489, 22509, 22529, 22549, 22569, 22589, 22609, 22629, 22649, 22669, 22689, 22709, 22729, 22749, 22769, 22789, 22809, 22829, 22849, 22869, 22889, 22909, 22929, 22949, 22969, 22989, 23017, 23037, 23057, 23077, 23097, 23117, 23137, 23157, 23177, 23197, 23217, 23237, 23257, 23277, 23297, 23317, 23337, 23364, 23385, 23405]
                xData = [6355.8, 6352.2, 6348.45, 6344.1, 6339.6, 6334.35, 6330.15, 6326.55, 6322.5, 6319.65, 6316.2, 6312.6, 6309, 6305.85, 6302.25, 6299.25, 6295.05, 6292.05, 6289.05, 6285.6, 6282, 6278.7, 6275.7, 6273.75, 6271.2, 6268.35, 6267.45, 6265.95, 6263.7, 6261.3, 6259.2, 6256.65, 6255.45, 6252.75, 6250.8, 6249.45, 6248.25, 6246, 6244.65, 6242.85, 6240.75, 6239.85, 6239.1, 6237.9, 6236.55, 6235.5, 6233.7, 6233.1, 6231.15, 6229.65, 6227.1, 6226.65, 6225, 6225.75, 6225.15, 6225.15, 6224.7, 6225.6, 6225.6, 6226.35, 6225.15, 6225.3, 6224.55, 6224.55, 6223.95, 6222.45, 6222.6, 6222.9, 6223.05, 6223.05, 6223.2, 6224.25, 6226.05, 6226.8, 6227.25, 6226.8, 6228.15, 6228.9, 6229.2, 6228.6, 6229.05, 6228.9, 6229.95, 6230.7, 6231.9, 6233.4, 6234.45, 6234.45, 6236.7, 6238.2, 6239.1, 6240, 6240.9, 6243.15, 6245.1, 6246.3, 6248.1, 6248.7, 6249.9, 6250.65, 6251.1, 6252.6, 6253.8, 6255.6, 6258.3, 6259.95, 6261, 6262.95, 6264.6, 6266.55, 6267.6, 6269.4, 6271.2, 6273.3, 6275.25, 6277.65, 6279.3, 6280.95, 6282, 6284.25, 6286.95, 6289.35, 6291.45, 6293.85, 6296.55, 6299.85, 6301.8, 6303.6, 6305.4, 6307.35, 6309.9, 6312.3, 6315, 6317.25, 6320.4, 6323.7, 6326.55, 6329.1, 6332.25, 6334.5, 6337.5, 6339.75, 6341.85, 6344.55, 6346.5, 6349.65, 6352.35, 6355.8, 6359.55, 6363, 6366, 6368.55, 6370.65, 6372.9, 6375, 6376.5, 6379.35, 6382.2, 6386.25, 6389.85, 6394.35, 6397.8, 6402.15, 6405.45, 6408, 6411.9, 6414.45, 6417.3, 6421.05, 6424.65, 6427.95, 6432.3, 6434.55, 6437.4, 6439.35, 6442.2, 6445.35, 6448.05, 6450, 6451.5, 6453.75, 6456.3, 6457.8, 6458.1, 6459.15, 6460.05, 6463.05, 6463.95, 6466.2, 6467.7, 6469.35, 6470.7, 6472.65, 6473.85, 6475.5, 6475.8, 6476.85, 6478.2, 6478.65, 6479.55, 6480, 6481.8, 6483, 6483, 6483.6, 6485.25, 6485.25, 6484.95, 6485.4, 6485.25, 6486.15, 6486.6, 6486, 6486.75, 6487.35, 6487.05, 6487.35, 6487.35, 6486.75, 6487.05, 6486.9, 6486.45, 6485.55, 6485.85, 6485.1, 6484.35, 6483.15, 6482.7, 6482.4, 6481.65, 6479.55, 6479.1, 6478.95, 6478.35, 6477, 6475.95, 6474.9, 6474.15, 6472.05, 6470.25, 6468.3, 6465.45, 6462.9, 6460.8, 6458.55, 6457.05, 6455.55, 6455.25, 6454.2, 6453.6, 6452.1, 6450.6, 6448.65, 6445.95, 6444.15, 6442.2, 6439.65, 6437.25, 6435.3, 6433.2, 6431.4, 6428.55, 6426.6, 6424.8, 6424.05, 6422.25, 6420.45, 6419.25, 6417.45, 6415.65, 6412.8, 6409.8, 6406.05, 6403.5, 6400.2, 6396, 6393.15, 6390.6, 6387.3, 6385.35, 6382.35, 6379.5, 6377.4, 6374.85, 6372.75, 6370.35, 6368.1, 6364.8, 6361.8, 6358.5, 6354.75, 6351.6, 6347.85, 6343.65, 6341.1, 6338.1, 6334.8, 6331.65, 6327.45, 6323.85, 6320.7, 6317.85, 6314.55, 6312, 6308.85, 6306.15, 6303.3, 6299.85, 6295.8, 6292.65, 6289.35, 6286.05, 6283.2, 6280.8, 6278.85, 6276.6, 6274.65, 6272.55, 6271.95, 6270.6, 6268.2, 6267.15, 6265.05, 6262.8, 6261.15, 6258.6, 6256.65, 6255.45, 6253.05, 6251.25, 6249.15, 6246.75, 6245.1, 6243.45, 6241.65, 6239.4, 6239.4, 6238.8, 6237.45, 6236.1, 6234.6]
                yData = [10752.45, 10752, 10751.25, 10750.8, 10751.25, 10750.35, 10750.5, 10749.9, 10751.1, 10752.15, 10752.9, 10753.8, 10754.7, 10754.85, 10755.75, 10755.6, 10755.45, 10754.85, 10754.7, 10755, 10755.3, 10755.3, 10755.15, 10755.75, 10756.35, 10756.8, 10757.85, 10758.6, 10759.95, 10761.3, 10762.35, 10763.25, 10764.75, 10765.05, 10767.15, 10767.75, 10769.25, 10770.75, 10773.3, 10773.9, 10776, 10776.6, 10778.25, 10779.15, 10779.9, 10780.05, 10782.3, 10783.5, 10785.3, 10786.5, 10787.85, 10788.3, 10789.5, 10790.55, 10791.9, 10793.1, 10793.4, 10794.6, 10796.55, 10797.3, 10798.05, 10798.35, 10798.95, 10800.75, 10802.4, 10803.6, 10805.25, 10806.45, 10808.7, 10809.6, 10811.25, 10812.45, 10813.8, 10815, 10815.75, 10816.35, 10818.15, 10818.45, 10818.45, 10819.8, 10821.3, 10822.35, 10823.25, 10824.3, 10825.5, 10826.55, 10827, 10827.75, 10829.55, 10830.45, 10831.2, 10832.4, 10834.5, 10836.15, 10836.9, 10838.1, 10840.2, 10841.1, 10842.15, 10842.75, 10844.25, 10845.3, 10845.45, 10846.05, 10846.65, 10846.8, 10848.75, 10849.2, 10850.4, 10851.6, 10852.5, 10854.3, 10855.8, 10855.05, 10856.4, 10857, 10857.3, 10857.6, 10858.8, 10859.55, 10860.15, 10859.1, 10859.55, 10860.75, 10861.95, 10861.8, 10862.4, 10863.15, 10864.2, 10865.4, 10865.25, 10866.3, 10867.65, 10868.4, 10869.75, 10870.5, 10870.05, 10871.4, 10871.1, 10870.2, 10870.65, 10870.05, 10870.65, 10871.55, 10871.85, 10872.75, 10873.35, 10873.5, 10873.8, 10874.4, 10874.85, 10874.25, 10873.5, 10873.8, 10873.8, 10873.65, 10873.2, 10872.9, 10873.05, 10872.6, 10873.05, 10872.75, 10873.05, 10872.45, 10872, 10871.25, 10871.55, 10870.5, 10869.3, 10868.25, 10867.8, 10866.75, 10867.05, 10865.1, 10864.2, 10863.6, 10863.3, 10861.95, 10861.2, 10860, 10861.2, 10860.15, 10859.55, 10858.95, 10858.5, 10857.45, 10856.85, 10854.3, 10854, 10852.65, 10850.25, 10848.3, 10847.4, 10845.9, 10844.4, 10842.15, 10841.55, 10841.4, 10841.25, 10839.3, 10837.5, 10837.5, 10836.45, 10834.95, 10833.75, 10831.5, 10830.75, 10830, 10828.2, 10828.05, 10826.85, 10824.9, 10824.6, 10823.7, 10822.95, 10821.15, 10819.2, 10818, 10816.95, 10815.15, 10814.25, 10812, 10810.8, 10807.8, 10805.7, 10804.5, 10803.3, 10802.1, 10801.35, 10800.75, 10800.75, 10800.45, 10798.8, 10797.45, 10795.35, 10793.25, 10791.6, 10790.25, 10789.35, 10787.25, 10785.75, 10784.55, 10783.8, 10782.6, 10781.55, 10780.35, 10779.15, 10777.95, 10777.5, 10775.55, 10774.5, 10773.3, 10770.9, 10770.3, 10769.85, 10768.2, 10767.9, 10766.85, 10765.35, 10764.3, 10764, 10763.7, 10763.25, 10762.5, 10762.35, 10761, 10760.7, 10759.65, 10757.7, 10756.2, 10755.15, 10753.2, 10752.6, 10752, 10751.55, 10751.1, 10750.5, 10749.75, 10749.45, 10750.65, 10750.35, 10750.65, 10750.05, 10749.9, 10749.6, 10749, 10747.95, 10747.05, 10746, 10746.15, 10744.8, 10744.8, 10744.35, 10743.75, 10744.5, 10744.65, 10744.5, 10744.5, 10744.5, 10745.55, 10745.25, 10744.2, 10744.65, 10744.8, 10745.85, 10746.3, 10746.15, 10747.35, 10747.8, 10748.1, 10748.4, 10749, 10749.3, 10749.45, 10750.05, 10751.7, 10751.25, 10752.45, 10753.2, 10752.9, 10753.8, 10754.25, 10754.7, 10755.6, 10755.75, 10757.25, 10759.05, 10759.95, 10761.15, 10761.45, 10763.55, 10764.9, 10764.9, 10766.1, 10767, 10767.9, 10770, 10770.45, 10771.35, 10773.6, 10775.1]
                zData = [2927.85, 2928.6, 2929.65, 2929.8, 2930.25, 2929.5, 2929.05, 2929.05, 2928.45, 2928.6, 2929.65, 2929.2, 2929.5, 2930.7, 2931.9, 2932.65, 2932.95, 2932.2, 2934.15, 2935.35, 2935.8, 2935.8, 2936.85, 2937.9, 2939.4, 2939.25, 2940.15, 2940, 2941.05, 2941.95, 2942.7, 2944.5, 2944.5, 2945.55, 2948.25, 2948.85, 2949.3, 2950.65, 2950.95, 2952.75, 2953.35, 2954.25, 2956.05, 2957.4, 2958.75, 2961.45, 2963.25, 2965.8, 2966.7, 2967.75, 2969.55, 2970.15, 2970, 2971.2, 2971.35, 2971.8, 2972.85, 2973.45, 2975.7, 2977.2, 2978.4, 2979.9, 2981.7, 2983.05, 2984.55, 2984.55, 2985.15, 2985.75, 2986.5, 2987.55, 2987.55, 2988.6, 2989.65, 2991.3, 2991.75, 2992.35, 2993.85, 2996.25, 2997, 2997.9, 2999.1, 3000.3, 3001.95, 3001.95, 3002.25, 3002.85, 3005.1, 3005.25, 3006.9, 3007.5, 3009.75, 3009.6, 3011.85, 3011.85, 3013.65, 3014.4, 3015, 3015, 3017.4, 3016.35, 3017.7, 3017.85, 3018.75, 3020.55, 3021.75, 3022.65, 3024.15, 3025.05, 3025.95, 3027.15, 3027.6, 3027.9, 3029.25, 3030.45, 3030.9, 3032.55, 3033.6, 3034.35, 3035.55, 3035.4, 3037.05, 3038.25, 3038.25, 3037.8, 3038.1, 3038.4, 3039.45, 3039.75, 3039.6, 3040.05, 3041.25, 3042.45, 3043.5, 3043.8, 3043.8, 3043.65, 3044.55, 3044.7, 3044.4, 3044.4, 3045.45, 3046.05, 3047.4, 3046.95, 3048.45, 3049.05, 3050.1, 3048.9, 3049.2, 3048.9, 3050.1, 3048.9, 3048.9, 3048, 3049.5, 3049.05, 3049.2, 3048.75, 3048.9, 3048.9, 3048.6, 3047.4, 3047.1, 3047.25, 3046.5, 3046.5, 3045.9, 3045.75, 3045.45, 3045, 3043.65, 3043.8, 3042.3, 3041.7, 3041.25, 3039.75, 3038.85, 3037.5, 3035.85, 3034.8, 3033.6, 3031.65, 3031.05, 3030, 3029.55, 3028.8, 3028.2, 3027.75, 3027.75, 3026.85, 3026.55, 3025.5, 3024.45, 3024.15, 3022.35, 3021, 3020.1, 3018.15, 3017.55, 3016.8, 3014.7, 3013.65, 3013.35, 3011.85, 3011.4, 3009.9, 3007.65, 3006.75, 3005.7, 3003, 3002.25, 3000, 2999.4, 2998.5, 2997.6, 2996.55, 2996.55, 2995.2, 2994.45, 2992.35, 2990.25, 2988.75, 2986.95, 2985.45, 2983.95, 2982, 2981.1, 2980.95, 2980.05, 2978.55, 2976.9, 2976.45, 2976.3, 2975.25, 2973.15, 2970.6, 2970.15, 2968.05, 2964.75, 2962.95, 2961.3, 2960.7, 2960.4, 2959.8, 2959.2, 2959.05, 2957.55, 2956.5, 2955.75, 2954.25, 2951.4, 2950.8, 2950.05, 2949.75, 2948.25, 2947.35, 2947.2, 2947.5, 2946.15, 2945.7, 2944.35, 2944.2, 2943.3, 2942.25, 2941.65, 2941.35, 2940.9, 2941.05, 2939.7, 2938.35, 2937.6, 2936.7, 2935.5, 2934.9, 2934.15, 2933.85, 2933.4, 2933.1, 2932.2, 2932.8, 2931.3, 2930.85, 2931.15, 2931.3, 2930.1, 2929.8, 2928.75, 2929.2, 2929.2, 2928.3, 2927.4, 2928.45, 2929.35, 2929.05, 2929.05, 2929.5, 2930.55, 2931.75, 2931.15, 2930.55, 2931.15, 2931.45, 2931, 2931, 2931.45, 2930.85, 2931.6, 2932.65, 2932.8, 2932.95, 2933.55, 2933.85, 2935.65, 2936.1, 2935.65, 2936.4, 2937.3, 2937.3, 2937.45, 2937.45, 2938.05, 2938.5, 2937.9, 2938.5, 2939.4, 2940, 2941.35, 2942.4, 2944.05, 2946.3, 2947.35, 2948.1, 2949.45, 2950.5, 2951.4, 2952.6, 2953.35, 2954.25, 2955.3, 2956.35, 2958.15]
                xTest = [6473.1, 6471.75, 6473.55, 6407.25, 6302.1, 6230.55, 6221.55, 6287.4, 6388.35]
                yTest = [10843.5, 10841.25, 10796.7, 10753.35, 10740, 10765.95, 10814.55, 10855.05, 10861.65]
                zTest = [3019.8, 3023.7, 2975.4, 2936.4, 2930.25, 2953.2, 2995.35, 3037.8, 3048]
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
        // do the same for the test cases
        for (let n = 0; n < xTest.length; n++) {
            let xyz = []
            xyz.push(xTest[n])
            xyz.push(yTest[n])
            xyz.push(zTest[n])
            testData.push(xyz)
        }
    }
}
