
 // TODO? No use is yet made of the accelerometer. Although seemingly helpful to compensate for
 // static tilt, on a moving buggy dynamic sideways accelerations confound the measurement of "down",
 // so applying tilt-compensation could actually make compass-heading accuracy worse!


// OPERATIONAL MODES (controlling data-logging)
enum Mode {
    Normal, // Normal usage, mounted in a buggy
    Capture, // Acquire a new test dataset, using a special rotating jig
    Analyse, // Test & debug (NOTE: named sets of test-dataset are hard-coded below)
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
    const EnoughScanTime = 1500 // minimum acceptable scan-time
    const EnoughSamples = 70 // fewest acceptable scan samples
    const TooManySamples = 500 // don't be too greedy with memory!
    const MarginalField = 10 // minimum acceptable field-strength for magnetometer readings
    const Circular = 1.03 // maximum eccentricity to consider an Ellipse as "circular" (3% gives ~1 degree error)
    const NearEnough = 0.75 // major-axis candidates must be longer than (longest*NearEnough)
                            // minor-axis candidates must be shorter than (shortest/NearEnough) 
    const Window = 7 // number of magnetometer samples needed to form a good average
    const SampleGap = 15 // minimum ms to wait between magnetometer readings
    const AverageGap = 25 // (achieved in practice, due to system interrupts)
    const Latency = Window * AverageGap // resulting time taken to collect a good moving average from scratch

    // SUPPORTING CLASSES

    /* =================================================================
    // An Arrow is an object holding a directed vector {u,v} in both Cartesian and Polar coordinates. 
    // It also carries a time-field, used to timestamp scanned samples.
    // An Arrow is used to hold a 2D magnetometer measurement as a re-centred vector.
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
    */
  
    // A Smoother object computes a moving average from a sequence of time-stamped values: 
    // in our case, magnetometer readings and their derivatives.
    // Timing irregularites due to scheduler interrupts demand this somewhat complex maths.
    // The constant {Window} governs the latency of the exponential averaging process.
    // Smoothers can work with arbitrary-sized vectors of values that share the same timestamp.
    // history[], previous[], latest[] and result[] arrays will be either 3-D, (for initial [X,Y,Z] scanning)
    // or 2-D, (for analysing chosen the [uDim,vDim]).
    class Smoother {
        dims: number; // dimensionality
        average: number[] = []; // 
        lastTime: number;
        lastInputs: number[] = [];

        constructor(first: number[], start: number) {
            this.dims = first.length
            this.lastTime = start
            for (let i = 0; i < this.dims; i++) {
                this.average.push(first[i])
                this.lastInputs.push(first[i])
            }
        }

        update(values: number[], timeStamp: number): number[] {
            // work out appropriate blend, based on time-step
            let timeFraction = (timeStamp - this.lastTime) / Latency
            let keepOld = Math.exp(-timeFraction)
            let inherited = (1 - keepOld) / timeFraction
            // we amplify the most recent sample's contribution to the inherited average
            let boostLast = (inherited - keepOld)
            let addNew = (1 - inherited)
            // (blending proportions keepOld + boostLast + addNew will always add up to 100%)
            // apply blending to all elements of old and new data arrays
            let result: number[] = []
            for (let i = 0; i < this.dims; i++) {
                result.push(keepOld * this.average[i]
                    + boostLast * this.lastInputs[i]
                    + addNew * values[i])
            }
            // update history for next time around
            this.lastTime = timeStamp
            this.average = result
            this.lastInputs = values

            return result
        }
    }
   
    /*
    An Ellipse is an object holding the characteristics of the (typically) elliptical
    view formed when projecting the scan Spin-Circle onto a 2-axis View-plane {XY, YZ or ZX} 
    This foreshortened view means that evenly spaced headings will appear bunched around the ends 
    of the ellipse. In the extreme "side-on" case, the ellipse collapses to just a straight line!
    
    To correct for foreshortening we must stretch the ellipse back into a circle, by rescaling
    parallel to the minor axis of the ellipse until it matches the major axis. 
    The ratio of major-axis to minor-axis gives the necessary scaling factor {eccentricity}.

    Depending on the exact mounting orientation of the microbit in the buggy, this ellipse may be
    tilted with respect to a particular view's axes, so correction then becomes a three stage process:
    1) rotate a new point [u,v] by {-tilt} so the minor-axis lines up with the V-axis.
    2) stretch the V coordinate by {eccentricity}
    3) rotate back by {tilt} to give the corrected [u,v] from which the heading can be derived.
  
    This correction can in theory be applied in any of the three views (unless exactly side on to
    the Spin-Circle), but the best accuracy is obtained from the most circular of the three views.
    Readings on a near-circular Ellipse are barely fore-shortened at all, so we can skip correction!

    So for each view the two important Ellipse properties we must derive are its {tilt} and {eccentricity}, 
    which first requires detection of its major and minor axes. The maths for fitting an ellipse to noisy
    2D data is both complex and fairly inaccurate. Luckily we make use of the orthogonal third dimension
    (the Normal) to give a simpler, faster solution.

    The three magnetometer readings are related by the formula:  {x^2 + y^2 + z^2 = B} (where B is the 
    constant magnetic field). So, for example in the XY plane, the ellipse radius {x^2 + y^2} is at a maximum
    (i.e. passing its major-axis) when the field aligns with this plane and the z-value is basically zero.
    Conversely the radius is at a minimum (its minor-axis) where the field points farthest from the plane,
    and the z-value changes from growing to shrinking (either poitive or negative). 
    
    The same holds true for the other two planes: the Normal helps us find the two axes. For each of these
    three mutually-orthogonal views, we re-label its coordinates as {u,v}, with {w} as the
    third (orthogonal) coordinate.

        // While scanning clockwise, the projection of the field direction will appear to trace out an
        // elliptical view of the anti-clockwise Spin-Circle. If viewed from "above" the polar angles of
        // successive readings (in radians) INCREASE (anti-clockwise); if viewed from "below" they DECREASE.


    */
    class Ellipse {
        plane: string // View name (just for debug)

        uDim: number // horizontal axis of this View
        vDim: number // vertical axis of this View
        uOff: number // horizontal offset needed to re-centre this Ellipse along the U-axis
        vOff: number // vertical offset needed to re-centre this Ellipse along the V-axis

        uHi: number // major-axis u contributions
        vHi: number // major=axis v contributions
        nHi: number // major-axis contributors
        rHi: number // major-axis magnitude (radius)

        uLo: number // minor-axis u contributions
        vLo: number // minor=axis v contributions
        nLo: number // minor-axis contributors
        rLo: number // minor-axis magnitude (radius)

        above: boolean  // flag saying {w} is currently "above" our plane
        gotMinor: boolean // flag saying we have clocked next transit of the minor-axis
        turns: number // number of full rotations since the first major-axis transit
        start: number // timestamp of first major-axis transit
        finish: number // timestamp of latest major-axis transit

        // calibration characteristics
        period: number // scan-rotation time as viewed in this plane
        tilt: number // direction of the Ellipse major-axis 
        cosa: number // helpful cosine of the tilt
        sina: number // helpful sine of the tilt
        eccentricity: number // ratio of major-axis to minor-axis magnitudes for this Ellipse
        isCircular: boolean // flag saying this "Ellipse" View is almost circular, simplifying future handling
        rotationSense: number // rotation sign = +/-1, reflecting this Ellipse's view of the clockwise scan


        constructor(plane: string, uDim: number, vDim: number, uOff: number, vOff:number) {
            this.plane = plane // just for debug...
            this.uDim = uDim
            this.vDim = vDim
            this.uOff = uOff
            this.vOff = vOff
            this.uHi = 0
            this.vHi = 0
            this.nHi = 0
            this.uLo = 0
            this.vLo = 0
            this.nLo = 0
            this.above = true
            this.gotMinor = false
            this.tilt = 0
            this.start = 0
            this.finish = 0
            this.turns = -1
            this.period = -1
            this.rotationSense = 1
        }

        // addMajor() is called whenever the Normal coordinate changes sign, indicating
        // that the Field is just crossing our plane and we're passing the major-axis.
        addMajor(i: number, u: number, v: number, w: number, t: number) {
            if (w > 0) { // just surfaced above our plane so add-in current vector coordinates
                this.above = true
                this.uHi += u
                this.vHi += v
                if (this.start == 0) this.start = t // start measuring turns
                this.finish = t
                this.turns++
            } else {  // just dipped below our plane so subtract current vector coordinates
                    // (as it's the other end of the major-axis)
                this.above = false
                this.uHi -= u
                this.vHi -= v
            }
            this.nHi++
            this.gotMinor = false

            if ((mode == Mode.Trace) || (mode == Mode.Analyse)) {
                datalogger.log(
                    datalogger.createCV("view", this.plane),
                    datalogger.createCV("index", i),
                    datalogger.createCV("uHi", round2(this.uHi)),
                    datalogger.createCV("vHi", round2(this.vHi)),
                    datalogger.createCV("nHi", round2(this.nHi)),
                    datalogger.createCV("uLo", round2(this.uLo)),
                    datalogger.createCV("vLo", round2(this.vLo)),
                    datalogger.createCV("nLo", round2(this.nLo)),
                    datalogger.createCV("above", this.above),
                    datalogger.createCV("gotMinor", this.gotMinor)
                )
            }
        }

        // addMinor() is called whenever the Normal coordinate changes from growing to shrinking.
        // This means the field-vector is farthest from our plane (above or below), so we're near
        // the minor-axis of this plane's ellipse.
        /* 
        NOTE: noisy readings may yield multiple inflections in this third coordinate 
         (e.g. peak-trough-peak above the plane, or trough-peak-trough below), leading to this method 
         being invoked multiple times for the same transit (though always an odd number).
         At peaks the radius is always added to the sum-vector [uLo, vLo] , while at troughs it is 
         always subtracted, so redundant matched pairs will always cancel out. However, we must be 
         careful only to count this transit once (using the flag this.newMinor), however many false 
         pairs of calls we may get.
        */
        addMinor(i: number, u: number, v: number, dw: number) {
            if (this.start > 0 ) { // don't start adding minor-axes until this.above is clearly known
                if (dw < 0) {  
                    this.uLo += u
                    this.vLo += v
                } else {
                    this.uLo -= u
                    this.vLo -= v
                }
                // because of the possibility of clustered peaks/troughs (some of which partially cancel) 
                // we only clock the first one on each transit
                if (!this.gotMinor) {
                    this.gotMinor = true
                    this.nLo++
                }
            }
            if ((mode == Mode.Trace) || (mode == Mode.Analyse)) {
                datalogger.log(
                    datalogger.createCV("view", this.plane),
                    datalogger.createCV("index", i),
                    datalogger.createCV("uHi", round2(this.uHi)),
                    datalogger.createCV("vHi", round2(this.vHi)),
                    datalogger.createCV("nHi", round2(this.nHi)),
                    datalogger.createCV("uLo", round2(this.uLo)),
                    datalogger.createCV("vLo", round2(this.vLo)),
                    datalogger.createCV("nLo", round2(this.nLo)),
                    datalogger.createCV("above", this.above),
                    datalogger.createCV("newMinor", this.gotMinor)
                )
            }
        }

        // calculate() method is called once all scandata has been processed
        calculate() {
            this.eccentricity = -1
            this.period = -1
            this.tilt = 0
            if (this.nHi > 0) {
                this.uHi /= this.nHi
                this.vHi /= this.nHi
                this.rHi = Math.sqrt(this.uHi * this.uHi + this.vHi * this.vHi)

                if (this.nLo > 0) {
                    this.uLo /= this.nLo
                    this.vLo /= this.nLo
                    this.rLo = Math.sqrt(this.uLo * this.uLo + this.vLo * this.vLo)
                    this.eccentricity = this.rHi / this.rLo

                    /* use vector addition to refine the axis tilt (but only
                    // after turning the minor-axis through a right angle)
                    let uMean = this.uHi + this.vLo
                    let vMean = this.vHi - this.uLo
                    let rMean = this.rHi + this.rLo
                    this.tilt = Math.atan2(vMean, uMean)
                    this.cosa = uMean / rMean
                    this.sina = vMean / rMean
                    */

                    // save major-axis and its components
                    this.tilt = Math.atan2(this.vHi, this.uHi)
                    this.cosa = this.uHi / this.rHi
                    this.sina = this.vHi / this.rHi
                }
            }
            if (this.turns > 0) {
                this.period = (this.finish - this.start) / this.turns
            }
        }
    }




    // GLOBALS

    let scanTimes: number[] = [] // sequence of time-stamps for scanned readings 
    let scanData: number[][] = [] // scanned sequence of [X,Y,Z] magnetometer readings
    let scanTime: number = 0 // duration of scan in ms
    //let views: Ellipse[] = [] // the three possible elliptical views of the Spin-Circle
    let bestView = -1
    let uDim = -1 // the "horizontal" axis (called U) for the best View
    let vDim = -1 // the "vertical" axis (called V) for the best View
    let north = 0 // reading registered as "North"
    let strength = 0 // the average magnetic field-strength observed by the magnetometer
    let period = -1 // overall assessment of average rotation time

    let xOff: number
    let yOff: number
    let zOff: number

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

    // new maths
    let xy: Ellipse
    let yz: Ellipse
    let zx: Ellipse


    // EXPORTED USER INTERFACES   

    /** 
     * Assuming the buggy is currently spinning clockwise on the spot, capture a time-stamped
     * sequence of magnetometer readings (TODO: and when finished, process them to set up the compass.)
     *
     * @param ms scanning-time in millisecs (long enough for more than one full rotation)    
     */

    //% block="scan clockwise for (ms) $ms" 
    //% inlineInputMode=inline 
    //% ms.shadow="timePicker" 
    //% ms.defl=0 
    //% weight=90 

    export function scanClockwise(ms: number) {
        // Sample magnetometer readings periodically over the specified duration (generally a couple
        // of seconds), and append a new [X,Y,Z] triple to the scanData[] array.
        // A timestamp for each sample is also recorded in the scanTimes[] array.

        // NOTE: To smooth out jitter, each reading is always a moving average of several consecutive readings.
        // Because sample-times may be irregular (due to scheduled interrupts), a Smoother is used to provide
        // a timing-aware exponential moving-average. The sample-grouping and spacing are controlled 
        // respectively by the constants Window and SampleGap, which together determine the Latency.

        scanTimes = []
        scanData = []
        let index = 0

        if (mode != Mode.Normal) {
            datalogger.deleteLog()
            datalogger.includeTimestamp(FlashLogTimeStampFormat.Milliseconds)
        }

        if (mode == Mode.Analyse) {
            simulateScan(dataset)
            basic.pause(ms)
        } else { // use live magnetometer
            let timeWas: number
            let timeNow: number
            let fresh: number[] = []
            let updated: number[] = []

            basic.pause(200) // wait for motors to stabilise (after initial kick-start)
            // get initial reading
            let timeStamp = input.runningTime()
            fresh = [
                input.magneticForce(Dimension.X),
                input.magneticForce(Dimension.Y),
                input.magneticForce(Dimension.Z)]

            // use a Smoother to maintain a rolling average
            let scan = new Smoother(fresh, timeStamp)

            // after an initial settling period, continue cranking out updated moving averages... 
            let startTime = timeStamp + Latency
            let stopTime = timeStamp + ms

            // ...until we run out of time (or space!)
            while ((timeStamp < stopTime)
                && (scanTimes.length < TooManySamples)) {
                // After processing, sleep until it's time for next sample.
                // NOTE: here is where various system subprograms will get scheduled.
                // If they need more time than we've offered, our next sample will get delayed!
                // (This seems to incur extra delays of ~44 ms every 100ms, plus ~26ms every 400ms)

                timeWas = timeStamp // remember time of latest sample
                timeNow = input.runningTime()
                basic.pause((timeWas + SampleGap) - timeNow) // pause for remainder of SampleGap (if any!)
                timeStamp = input.runningTime() // take a fresh set of readings

                fresh = [
                    input.magneticForce(Dimension.X),
                    input.magneticForce(Dimension.Y),
                    input.magneticForce(Dimension.Z)]
                updated = scan.update(fresh, timeStamp)

                // only start recording once the moving average has stabilised
                if (timeStamp > startTime) {
                    // store the triple of averaged [X,Y,Z] values (as a deep copy!)
                    scanData.push([updated[0], updated[1], updated[2]])
                    scanTimes.push(timeStamp)  // timestamp it

                    if (mode != Mode.Normal) {
                        datalogger.log(
                            datalogger.createCV("t", timeStamp),
                            datalogger.createCV("x", round2(updated[0])),
                            datalogger.createCV("y", round2(updated[1])),
                            datalogger.createCV("z", round2(updated[2]))
                        )
                    }
                }
                index++
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
        // reset global defaults
        bestView = -1
        strength = -1
        period = -1

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
        // The means of the extremes give an approximation to the central offsets.
        xOff = (xhi + xlo) / 2
        yOff = (yhi + ylo) / 2
        zOff = (zhi + zlo) / 2

        // create three Ellipse instances, for analysing each possible 2D view of the spin-Circle
        xy = new Ellipse("XY", Dimension.X, Dimension.Y, xOff, yOff)
        yz = new Ellipse("YZ", Dimension.Y, Dimension.Z, yOff, zOff)
        zx = new Ellipse("ZX", Dimension.Z, Dimension.X, zOff, xOff)

        extractAxes()

        // check that at least one View saw at least one complete rotation (with a measurable period)...
        if ((xy.period + yz.period + zx.period) < 0 ) {
            return -3 // "NOT ENOUGH SCAN ROTATION"
        }

        // Choose the "roundest" Ellipse view  --the one with lowest eccentricity.
        let view: Ellipse = xy
        if (yz.eccentricity < view.eccentricity) view = yz
        if (zx.eccentricity < view.eccentricity) view = zx

        // periodicity is unreliable in a near-circular View: average just the other two views' measurements
        period = (xy.period + yz.period + zx.period - view.period) / 2

        // For efficiency, extract various characteristics from our adopted best Ellipse
        uDim = view.uDim
        vDim = view.vDim
        uOff = view.uOff
        vOff = view.vOff
        scale = view.eccentricity // the scaling needed to balance axes
        isCircular = (scale <= Circular)
        theta = view.tilt // the rotation (in radians) of the major-axis from the U-axis
        cosTheta = view.cosa
        sinTheta = view.sina
        
        // Having successfully set up the projection parameters for the bestView, get a
        // stable fix on the current heading, which we will then designate as "North".
        // (This is the global fixed bias to be subtracted from all future readings)
        north = 0
        north = takeSingleReading()

        if ((mode == Mode.Trace) || (mode == Mode.Analyse || (mode == Mode.Capture))) {
            datalogger.log(
                datalogger.createCV("view", view.plane),
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
        //views = []

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
    // Depending on mounting orientation, the bestView might possibly be seeing the Spin-Circle from
    // "underneath", with the field-vector appearing to move clockwise  --effectively experiencing an
    // anti-clockwise scan. In this case the rotationSense will be negative.
        let heading = (takeSingleReading() - north) * rotationSense
        // NOTE: that there is a double reversal going on here:
        // Viewed from above, the Field-vector reading in radians INCREASES (anticlockwise) w.r.t "North"
        // as the buggy's compass-heading INCREASES (clockwise).
        // The reading will therefore increase by HalfPi after a right-turn from initially facing "North".
        // So subtracting North (cyclically), that converts asDegrees() to +90 
        // From below (when rotationSense = -1), the same right-turn would DECREASE the reading by HalfPi
        // necessitating a third reversal, but after first subtracting North !
        return asDegrees(heading)
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
        if (period < 0) {
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
        //views = []
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

    /** Take the sum of several new readings to get a stable fix on the current heading.
     *  @return the angle of the magnetic field-vector (in radians anticlockwise
     * from the horizontal U-axis), corrected for any fore-shortening due to projection
     * onto the bestView plane.
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
        switch (mode) {
            case Mode.Normal:
                // get a new sample as the average of {Window} consecutive 2D readings, {SampleGap} apart
                for (let i = 0; i < Window; i++) {
                    basic.pause(SampleGap)
                    uRaw += input.magneticForce(uDim)
                    vRaw += input.magneticForce(vDim)
                }
                uRaw /= Window
                vRaw /= Window
                break

            case Mode.Analyse: // just choose the next test-data value
                uRaw = testData[test][uDim]
                vRaw = testData[test][vDim]
                test = (test + 1) % testData.length
                break

            case Mode.Capture: // capture averages of {Window} consecutive 3D readings, {SampleGap} apart
            case Mode.Trace:
                let xyz = [0, 0, 0]
                for (let i = 0; i < Window; i++) {
                    basic.pause(SampleGap)
                    xyz[Dimension.X] += input.magneticForce(Dimension.X)
                    xyz[Dimension.Y] += input.magneticForce(Dimension.Y)
                    xyz[Dimension.Z] += input.magneticForce(Dimension.Z)
                }
                xyz[Dimension.X] /= Window
                xyz[Dimension.Y] /= Window
                xyz[Dimension.Z] /= Window

                datalogger.setColumnTitles("t", "x", "y", "z")
                datalogger.log(
                    datalogger.createCV("t", input.runningTime()),
                    datalogger.createCV("x", round2(xyz[Dimension.X])),
                    datalogger.createCV("y", round2(xyz[Dimension.Y])),
                    datalogger.createCV("z", round2(xyz[Dimension.Z])))

                // now pick the coordinates we want for the current view
                uRaw = xyz[uDim]
                vRaw = xyz[vDim]
                break
        }

        // re-centre this latest point w.r.t our Ellipse origin
        u = uRaw - uOff
        v = vRaw - vOff

        if (isCircular) {
            reading = Math.atan2(v, u)
        } else {
            // Unless this Ellipse.isCircular, any {u,v} reading will be foreshortened in this View, and
            // must be stretched along the Ellipse minor-axis to place it correctly onto the Spin-Circle.

            // We must first rotate our point CLOCKWISE (by -theta) to align the Ellipse minor-axis 
            // angle with the V-axis 
            uNew = u * cosTheta + v * sinTheta
            vNew = v * cosTheta - u * sinTheta

            // Now scale up, just along V, re-balancing the axes to make the Ellipse circular
            // thereby moving our point outwards onto the rim of the Spin-Circle
            uFix = uNew
            vFix = vNew * scale
            // get the adjusted angle for this corrected {u,v}
            reading = Math.atan2(vFix, uFix)
            // finally, undo our first rotation by adding theta back in
            reading = (reading + theta) % TwoPi
        }

        if ((mode == Mode.Trace) || (mode == Mode.Analyse)) {
            // just for debug, show coordinates of "stretched" reading after undoing rotation
            let uRim = uFix * cosTheta - vFix * sinTheta
            let vRim = vFix * cosTheta + uFix * sinTheta
            datalogger.log(
                datalogger.createCV("uNew", round2(uNew)),
                datalogger.createCV("vNew", round2(vNew)),
                datalogger.createCV("uFix", round2(uFix)),
                datalogger.createCV("vFix", round2(vFix)),
                datalogger.createCV("uRim", round2(uRim)),
                datalogger.createCV("vRim", round2(vRim)),
                datalogger.createCV("reading", round2(reading)),
                datalogger.createCV("[reading]", round2(asDegrees(reading) * rotationSense)),
                datalogger.createCV("bearing", round2(reading - north)),
                datalogger.createCV("[bearing]", round2(asDegrees(reading - north) * rotationSense))
            )
        }
        return reading
    }

    /* Use vector addition to average a sequence of candidate axis Arrows (conceivably none!)
    // The ones pointing away from the evolving average (even slightly) are assumed to belong
    // to the other end of the axis, so will get reversed.
    // The returned Arrow shows the average axis length and angle.
    // Assuming candidates represent more than one revolution, the periodocity is also calculated.
    
    function computeAxis(sheaf: Arrow[]): Arrow {
        let result = new Arrow(0, 0, 0)
        let turns = 0
        let startTime = 0
        let endTime = 0
        let period = -1
        let flipped = false
        let uSum = 0
        let vSum = 0
        let rSum = 0
        let count = sheaf.length
        if (count > 0) {
            // initialise axis as first (or only?) candidate
            uSum = sheaf[0].u
            vSum = sheaf[0].v
            rSum = sheaf[0].size
            let axis = sheaf[0].angle
            startTime = sheaf[0].time
            flipped = false
            for (let i = 1; i < count; i++) {
                // does next candidate point nearer the head or the tail of the axis?
                if (Math.abs(angleSpan(axis, sheaf[i].angle)) < HalfPi) {
                    // chain this candidate onto the emerging axis
                    uSum += sheaf[i].u
                    vSum += sheaf[i].v
                    // the first unflipped candidate after one or more flipped ones clocks a new revolution
                    if (flipped) {
                        flipped = false
                        turns++
                        endTime = sheaf[i].time
                    }
                } else { // flip this arrow before chaining it, as it's pointing the "wrong" way
                    flipped = true
                    uSum -= sheaf[i].u
                    vSum -= sheaf[i].v
                }
                // get the new blended angle
                axis = Math.atan2(vSum, uSum)
                rSum += sheaf[i].size

                if ((mode == Mode.Trace) || (mode == Mode.Analyse)) {
                    datalogger.log(
                        datalogger.createCV("time", sheaf[i].time),
                        datalogger.createCV("uCand", round2(sheaf[i].u)),
                        datalogger.createCV("vCand", round2(sheaf[i].v)),
                        datalogger.createCV("rCand", round2(sheaf[i].size)),
                        datalogger.createCV("[aCand]", round2(asDegrees(sheaf[i].angle))),
                        datalogger.createCV("flip?", flipped),
                        datalogger.createCV("uSum", round2(uSum)),
                        datalogger.createCV("vSum", round2(vSum)),
                        datalogger.createCV("rSum", round2(rSum)),
                        datalogger.createCV("[axis]", round2(asDegrees(axis))))
                }
            }

            // build the result Arrow
            result.size =  rSum / count  // the average radius
            result.angle = axis
            result.u = result.size * Math.cos(axis)
            result.v = result.size * Math.sin(axis)

            // compute the average rotation time (so long as we've made at least one complete revolution)
            if (endTime > 0) {
                result.time = (endTime - startTime) / turns
            }
        }
        // hijack the time property to return the estimated period
        return result
    }
    */

    
    /** gives the signed difference between angles a & b (allowing for roll-round)
     * @param a first angle in radians
     * @param b second angle in radians
     * @returns the acute (i.e. smaller) difference in angle
     */
    function angleSpan(a:number, b:number){
        return((ThreePi + b - a) % TwoPi) - Math.PI
    }


    // helpful for logging...
    function round2(v: number): number {
        return (Math.round(100 * v) / 100)
    }


    // Convert an angle measured in radians to degrees.
    function asDegrees(angle: number): number {
        return ((angle * RadianDegrees) + 360) % 360
    }

    function extractAxes() {
        /*
        Re-centre all of the scanData samples, so eliminating "hard-iron" magnetic effects. 
        This will also re-centre the elliptical projections of the Spin-Circle in each 2D view.
        
        For correction of future readings we will need to find the eccentricity and tilt of each ellipse,
        and then (for highest accuracy) choose the most circular view.

        Detecting peaks and troughs in noisy data is error-prone, so we use a Smoother to minimise 
        (but not entirely eliminate) spurious inflection-points. Due to the latency of this moving average,
        (given by the constant Window) the minor-axis was actually passed {Window} samples earlier 
        than its actual detection.
        */

        let radiusUV = 0
        let radiusYZ = 0
        let radiusZX = 0
        let nXY = 0
        let nYZ = 0
        let nZX = 0

        let xWas = 0
        let yWas = 0
        let zWas = 0
        let x = scanData[0][Dimension.X] - xOff
        let y = scanData[0][Dimension.Y] - yOff
        let z = scanData[0][Dimension.Z] - zOff
        let t = scanTimes[0]
        let dx = 0
        let dy = 0
        let dz = 0
        let dxWas = 0
        let dyWas = 0
        let dzWas = 0
        let delta = new Smoother([dx, dy, dz], t)
        let delay: number[][] = []


        for (let i = 1; i < scanTimes.length; i++) {
            // update history
            xWas = x
            yWas = y
            zWas = z
            dxWas = dx
            dyWas = dy
            dzWas = dz

            // re-centre the next scanData sample
            x = scanData[i][Dimension.X] - xOff
            y = scanData[i][Dimension.Y] - yOff
            z = scanData[i][Dimension.Z] - zOff
            t = scanTimes[i]

            delay.push([x, y, z])  // this rolling array remembers recent sample history...


            // use the Smoother to get less noisy slopes
            delta.update([x - xWas, y - yWas, z - zWas], t)
            dx = delta.average[Dimension.X]
            dy = delta.average[Dimension.Y]
            dz = delta.average[Dimension.Z]

            // to aid detection of change of sign, we doctor any values that are exactly zero
            if (x == 0) x = xWas
            if (y == 0) y = yWas
            if (z == 0) z = zWas
            if (dx == 0) dx = dxWas
            if (dy == 0) dy = dyWas
            if (dz == 0) dz = dzWas

            // crossing a plane implies we're passing the major-axis of its ellipse
            if (z * zWas < 0) {
                xy.addMajor(i, x, y, z, t)
            }
            if (x * xWas < 0) {
                yz.addMajor(i, y, z, x, t)
            }
            if (y * yWas < 0) {
                zx.addMajor(i, z, x, y, t)
            }
            // don't start checking for minor axis until delta Smoother has stabilised
            if (i > Window) {
                delay.splice(0, 1) // always maintain exactly {Window} samples in the queue 
                // in case of minor-axis detections, recover readings for scanData[i-Window] 
                let xOld = delay[0][Dimension.X]
                let yOld = delay[0][Dimension.Y]
                let zOld = delay[0][Dimension.Z]
            
                if (dz * dzWas < 0) {
                    xy.addMinor(i, xOld, yOld, dz)
                }
                if (dx * dxWas < 0) {
                    yz.addMinor(i, yOld, zOld, dx)
                }
                if (dy * dyWas < 0) {
                    zx.addMinor(i, zOld, xOld, dy)
                }
            }
        }

        // use collected vector-sums to compute ellipse characteristics
        xy.calculate()
        yz.calculate()
        zx.calculate()

        // compare cross-products of last two samples to detect rotation-sense for each plane
        xy.rotationSense = (dx * dyWas) > (dy * dxWas) ? 1 : -1
        yz.rotationSense = (dy * dzWas) > (dz * dyWas) ? 1 : -1
        xy.rotationSense = (dz * dxWas) > (dx * dzWas) ? 1 : -1
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
            case "blup70": // angled with bottom-left upwards; dip = 70
                scanTimes = [31881, 31929, 31945, 31961, 31977, 31993, 32009, 32025, 32081, 32097, 32113, 32129, 32145, 32161, 32177, 32193, 32209, 32265, 32281, 32297, 32313, 32329, 32345, 32361, 32377, 32393, 32449, 32465, 32481, 32497, 32513, 32529, 32545, 32561, 32577, 32661, 32677, 32693, 32709, 32725, 32741, 32757, 32773, 32789, 32845, 32861, 32877, 32893, 32909, 32925, 32941, 32957, 32973, 33029, 33045, 33061, 33077, 33093, 33109, 33125, 33141, 33157, 33217, 33233, 33249, 33265, 33281, 33297, 33313, 33329, 33413, 33429, 33445, 33461, 33477, 33493, 33509, 33525, 33541, 33601, 33617, 33633, 33649, 33665, 33681, 33697, 33713, 33769, 33785, 33801, 33817, 33833, 33849, 33865, 33881, 33937, 33953, 33969, 33985, 34001, 34017, 34033, 34049, 34065, 34149, 34165, 34181, 34201, 34217, 34233, 34249, 34265, 34325, 34341, 34357, 34373, 34389, 34405, 34421, 34437, 34453, 34509, 34525, 34541, 34557, 34573, 34589, 34605, 34621, 34677, 34693, 34709, 34725, 34741, 34757, 34773, 34789, 34873, 34889, 34905, 34921, 34937, 34953, 34969, 34985, 35001, 35057, 35073, 35089, 35105, 35121, 35137, 35153, 35169, 35229, 35245, 35261, 35277, 35293, 35309, 35325, 35341, 35397, 35413, 35429, 35445, 35461, 35477, 35493, 35509, 35525, 35613, 35629, 35645, 35661, 35677, 35693, 35713, 35729, 35749, 35805, 35821, 35837, 35853, 35869, 35885, 35901, 35917, 35933, 35989, 36005, 36025, 36041, 36057, 36073, 36089, 36105, 36121, 36181, 36197, 36213, 36229, 36249, 36265, 36281, 36297, 36385, 36401, 36417, 36433, 36449, 36465, 36481, 36501, 36517, 36577, 36593, 36609, 36625, 36641, 36657, 36673, 36689, 36709, 36765, 36781, 36797, 36813, 36829, 36845, 36861, 36877, 36897, 36953, 36969, 36985, 37001, 37017, 37033, 37049, 37065, 37085, 37169, 37185, 37201, 37217, 37233, 37249, 37265, 37285, 37301, 37361, 37377, 37393, 37409, 37425, 37441, 37461, 37477, 37537, 37553, 37569, 37585, 37601, 37617, 37637, 37653, 37713]
                xData = [5.61, 5.99, 6.1, 6.24, 6.38, 6.48, 6.58, 6.65, 6.86, 6.91, 6.94, 6.97, 6.99, 7.02, 7.02, 6.97, 6.94, 6.76, 6.66, 6.6, 6.53, 6.46, 6.4, 6.33, 6.25, 6.16, 5.69, 5.55, 5.45, 5.31, 5.1, 4.92, 4.75, 4.57, 4.38, 3.22, 2.98, 2.76, 2.58, 2.43, 2.25, 2.02, 1.78, 1.52, 0.57, 0.3, 0.05, -0.2, -0.49, -0.84, -1.18, -1.48, -1.76, -2.85, -3.21, -3.62, -4.01, -4.36, -4.73, -5.08, -5.4, -5.78, -7.32, -7.69, -8.02, -8.39, -8.77, -9.14, -9.54, -9.99, -12.4, -12.8, -13.17, -13.51, -13.85, -14.22, -14.6, -14.94, -15.24, -16.55, -16.93, -17.3, -17.66, -18.03, -18.37, -18.67, -18.94, -20.03, -20.35, -20.66, -20.95, -21.22, -21.48, -21.75, -22.06, -23.01, -23.23, -23.44, -23.63, -23.82, -23.99, -24.15, -24.39, -24.6, -25.37, -25.48, -25.61, -25.77, -25.83, -25.87, -25.89, -25.94, -26.1, -26.12, -26.14, -26.14, -26.13, -26.08, -26, -25.94, -25.88, -25.66, -25.62, -25.53, -25.35, -25.15, -24.95, -24.75, -24.58, -23.88, -23.68, -23.47, -23.25, -23.07, -22.84, -22.58, -22.32, -20.76, -20.4, -20.01, -19.64, -19.27, -18.96, -18.64, -18.26, -17.86, -16.45, -16, -15.54, -15.12, -14.71, -14.31, -13.88, -13.45, -11.93, -11.58, -11.21, -10.83, -10.41, -10.02, -9.64, -9.25, -7.92, -7.55, -7.18, -6.81, -6.5, -6.21, -5.92, -5.59, -5.26, -3.36, -2.98, -2.61, -2.25, -1.88, -1.52, -1.07, -0.72, -0.3, 0.78, 1.08, 1.43, 1.79, 2.08, 2.32, 2.59, 2.87, 3.1, 3.81, 4.02, 4.31, 4.56, 4.81, 5.01, 5.18, 5.35, 5.5, 5.94, 6.03, 6.13, 6.19, 6.24, 6.28, 6.31, 6.31, 6.05, 5.99, 5.92, 5.85, 5.8, 5.75, 5.68, 5.56, 5.45, 5, 4.83, 4.64, 4.42, 4.16, 3.91, 3.72, 3.49, 3.22, 2.38, 2.11, 1.85, 1.57, 1.3, 1.03, 0.76, 0.51, 0.19, -0.83, -1.16, -1.46, -1.77, -2.08, -2.41, -2.71, -2.97, -3.32, -4.98, -5.31, -5.64, -5.95, -6.26, -6.52, -6.78, -7.18, -7.51, -8.79, -9.12, -9.42, -9.77, -10.13, -10.42, -10.8, -11.15, -12.33, -12.62, -12.9, -13.18, -13.47, -13.78, -14.14, -14.41, -15.47]
                yData = [48.69, 47.65, 47.33, 47.04, 46.68, 46.24, 45.78, 45.37, 43.94, 43.53, 43.13, 42.72, 42.33, 41.95, 41.54, 41.11, 40.67, 39.09, 38.66, 38.26, 37.85, 37.45, 37.04, 36.61, 36.2, 35.83, 34.61, 34.21, 33.79, 33.4, 33.04, 32.7, 32.36, 32.05, 31.74, 30, 29.65, 29.29, 28.99, 28.74, 28.48, 28.19, 27.86, 27.53, 26.5, 26.25, 25.97, 25.65, 25.37, 25.08, 24.76, 24.49, 24.21, 23.3, 23.09, 22.89, 22.68, 22.5, 22.33, 22.12, 21.91, 21.72, 21.05, 20.91, 20.8, 20.68, 20.54, 20.37, 20.22, 20.12, 19.75, 19.71, 19.68, 19.63, 19.61, 19.6, 19.59, 19.59, 19.61, 19.84, 19.94, 20.08, 20.24, 20.38, 20.47, 20.55, 20.67, 21.17, 21.33, 21.49, 21.64, 21.81, 22, 22.15, 22.34, 23.12, 23.34, 23.56, 23.75, 24, 24.3, 24.59, 24.85, 25.14, 26.87, 27.21, 27.55, 27.98, 28.32, 28.65, 28.99, 29.32, 30.62, 30.99, 31.37, 31.74, 32.11, 32.49, 32.85, 33.22, 33.63, 35.08, 35.49, 35.91, 36.36, 36.8, 37.26, 37.75, 38.18, 39.57, 39.97, 40.36, 40.77, 41.18, 41.55, 41.92, 42.34, 44.26, 44.63, 45, 45.32, 45.67, 46.03, 46.34, 46.68, 47.04, 48.06, 48.35, 48.65, 48.94, 49.16, 49.37, 49.61, 49.83, 50.48, 50.62, 50.75, 50.9, 51.08, 51.23, 51.37, 51.5, 51.8, 51.84, 51.94, 52.03, 52.05, 52.07, 52.07, 52.05, 52.03, 51.85, 51.83, 51.8, 51.7, 51.59, 51.47, 51.31, 51.19, 51.03, 50.39, 50.15, 49.88, 49.64, 49.41, 49.18, 48.96, 48.67, 48.33, 47.07, 46.68, 46.18, 45.75, 45.33, 44.95, 44.56, 44.16, 43.77, 42.35, 41.93, 41.47, 41.06, 40.61, 40.22, 39.78, 39.36, 37.26, 36.9, 36.52, 36.14, 35.77, 35.4, 35.01, 34.54, 34.18, 32.78, 32.41, 32.02, 31.57, 31.14, 30.75, 30.39, 30.04, 29.57, 28.43, 28.15, 27.83, 27.5, 27.18, 26.85, 26.53, 26.22, 25.85, 24.82, 24.55, 24.26, 23.99, 23.78, 23.53, 23.27, 23.08, 22.85, 21.83, 21.67, 21.54, 21.38, 21.22, 21.1, 20.98, 20.81, 20.67, 20.17, 20.1, 20.07, 19.98, 19.85, 19.71, 19.59, 19.51, 19.32, 19.3, 19.28, 19.28, 19.28, 19.26, 19.22, 19.24, 19.33]
                zData = [61.37, 61.87, 62.04, 62.23, 62.42, 62.63, 62.83, 63.06, 63.86, 64.03, 64.23, 64.42, 64.55, 64.71, 64.92, 65.12, 65.28, 65.85, 66.01, 66.19, 66.37, 66.5, 66.62, 66.79, 67, 67.17, 67.57, 67.6, 67.65, 67.76, 67.85, 67.98, 68.09, 68.18, 68.28, 68.64, 68.7, 68.69, 68.7, 68.76, 68.85, 68.93, 68.96, 69, 69.08, 69.06, 69.04, 69.07, 69.12, 69.1, 69.04, 68.99, 68.97, 68.93, 68.88, 68.83, 68.78, 68.72, 68.67, 68.59, 68.53, 68.45, 68.14, 68.05, 67.94, 67.83, 67.72, 67.64, 67.51, 67.37, 66.72, 66.56, 66.36, 66.15, 65.99, 65.86, 65.71, 65.54, 65.35, 64.65, 64.45, 64.26, 64.04, 63.83, 63.63, 63.42, 63.21, 62.46, 62.23, 62, 61.79, 61.57, 61.36, 61.17, 60.99, 60.25, 60.03, 59.82, 59.6, 59.4, 59.19, 58.95, 58.76, 58.58, 57.42, 57.2, 56.95, 56.69, 56.55, 56.37, 56.17, 56.02, 55.48, 55.28, 55.09, 54.95, 54.78, 54.57, 54.38, 54.24, 54.13, 53.76, 53.6, 53.44, 53.29, 53.14, 53.02, 52.87, 52.73, 52.28, 52.17, 52.06, 51.99, 51.96, 51.92, 51.87, 51.8, 51.47, 51.46, 51.49, 51.5, 51.49, 51.49, 51.52, 51.59, 51.67, 51.91, 51.98, 52.04, 52.07, 52.09, 52.15, 52.23, 52.3, 52.54, 52.61, 52.67, 52.75, 52.9, 53.11, 53.31, 53.45, 53.9, 54.02, 54.11, 54.21, 54.39, 54.59, 54.74, 54.9, 55.08, 56.06, 56.23, 56.4, 56.62, 56.83, 57.03, 57.32, 57.57, 57.79, 58.33, 58.53, 58.74, 58.97, 59.24, 59.53, 59.77, 59.95, 60.14, 61.02, 61.32, 61.66, 61.92, 62.22, 62.5, 62.76, 63, 63.26, 64.16, 64.32, 64.43, 64.61, 64.92, 65.13, 65.31, 65.49, 66.29, 66.43, 66.59, 66.76, 66.91, 67.05, 67.15, 67.26, 67.37, 67.75, 67.89, 68.03, 68.14, 68.26, 68.35, 68.4, 68.44, 68.51, 68.71, 68.75, 68.75, 68.77, 68.82, 68.87, 68.94, 69.02, 69.07, 69.11, 69.1, 69.08, 69.05, 69.02, 68.98, 68.98, 69, 68.97, 68.57, 68.5, 68.46, 68.38, 68.29, 68.23, 68.18, 68.13, 68.05, 67.69, 67.58, 67.45, 67.35, 67.25, 67.16, 67.06, 66.99, 66.64, 66.54, 66.44, 66.35, 66.29, 66.23, 66.1, 65.99, 65.53]
                xTest = [-25.08, -26.12, -17.66, -4.49, 4.99, 6.55, -2.46, -15.36, -25.47, -26.06, -18.31, -4.9, 5.21, 6.13, -2.34, -16.07, -24.82, -26.35, -18.15, -5.31, 4.85, 6.33, -1.89, -15.62, -25.71, -26.65, -18.33, -5, 4.72, 5.88, -2.5, -16.29, -25.92, -26.85, -18.47, -5.35, 4.68, 5.7, -2.68, -16.25, -25.69, -26.87, -18.45, -5.47, 4.68, 5.94, -2.52, -15.97, -25.31, -26.35, -18.49, -5.12, 5.05, 5.68, -2.68, -16.13, -25.8, -27.42, -18.41, -5.49, 4.76, 5.41, -2.83, -16.38, -26.04]
                yTest = [24.46, 36.27, 48.55, 52.71, 46.86, 34.78, 22.73, 18.1, 23.43, 35.59, 47.44, 52.37, 46.57, 34.56, 22.83, 18.65, 24.02, 35.91, 47.7, 52.06, 46.71, 33.93, 22.4, 17.92, 22.93, 35.22, 47.14, 51.82, 46.07, 33.57, 21.92, 17.31, 22.83, 34.17, 46.61, 51.12, 45.48, 33.55, 21.03, 17.11, 22.71, 34.46, 46.43, 51.28, 45.5, 33.17, 21.49, 16.97, 22.3, 33.79, 46.49, 51.32, 45.58, 32.84, 21.23, 16.45, 22.65, 34.32, 46.11, 50.39, 45.02, 32.68, 21.43, 16.73, 21.88]
                zTest = [58.77, 52.87, 51.14, 54.79, 61.46, 67.49, 69.05, 65.55, 58.29, 52.56, 50.62, 54.1, 61.38, 67.57, 69.12, 65.63, 58.54, 52.35, 50.72, 53.85, 60.86, 67.34, 68.74, 65.4, 58.67, 52.74, 51.22, 54.64, 61.25, 67.74, 69.49, 65.61, 58.98, 52.81, 50.95, 54.89, 61.4, 67.24, 69.39, 65.57, 58.54, 52.45, 51.43, 54.79, 61.34, 67.8, 69.47, 65.88, 58.4, 52.93, 51.7, 55.41, 61.38, 67.34, 69.8, 65.8, 58.75, 52.77, 50.99, 54.6, 61.69, 67.36, 69.39, 65.68, 58.42]
                break


            case "blup70_0714_1743":
                scanTimes = [32009, 32057, 32073, 32089, 32105, 32121, 32137, 32193, 32209, 32225, 32241, 32257, 32273, 32289, 32305, 32361, 32377, 32393, 32409, 32425, 32441, 32457, 32473, 32529, 32545, 32561, 32577, 32593, 32609, 32625, 32713, 32729, 32745, 32761, 32777, 32793, 32809, 32825, 32885, 32901, 32917, 32933, 32949, 32965, 32981, 33037, 33053, 33069, 33085, 33101, 33117, 33133, 33149, 33205, 33221, 33237, 33253, 33269, 33285, 33301, 33385, 33401, 33417, 33433, 33449, 33465, 33481, 33497, 33553, 33569, 33585, 33601, 33617, 33633, 33649, 33665, 33721, 33737, 33753, 33769, 33785, 33801, 33817, 33873, 33889, 33905, 33921, 33937, 33953, 33969, 33985, 34069, 34085, 34101, 34117, 34133, 34149, 34165, 34193, 34225, 34241, 34257, 34273, 34289, 34305, 34321, 34381, 34397, 34413, 34429, 34445, 34461, 34477, 34493, 34549, 34565, 34581, 34597, 34613, 34629, 34645, 34729, 34745, 34761, 34777, 34793, 34809, 34825, 34841, 34897, 34913, 34929, 34945, 34961, 34977, 34993, 35049, 35065, 35081, 35097, 35113, 35129, 35145, 35161, 35217, 35233, 35249, 35265, 35281, 35297, 35313, 35329, 35413, 35429, 35445, 35461, 35477, 35493, 35509, 35565, 35581, 35597, 35613, 35629, 35645, 35661, 35677, 35733, 35749, 35765, 35781, 35797, 35813, 35829, 35889, 35905, 35921, 35941, 35957, 35977, 35993, 36009, 36093, 36109, 36125, 36141, 36157, 36173, 36189, 36205, 36265, 36285, 36301, 36317, 36333, 36349, 36365, 36425, 36441, 36457, 36473, 36489, 36509, 36525, 36541, 36601, 36617, 36633, 36649, 36665, 36681, 36697, 36717, 36801, 36817, 36833, 36849, 36865, 36881, 36897, 36957, 36973, 36993, 37009, 37025, 37041, 37057, 37073, 37133, 37149, 37165, 37185, 37201, 37217, 37233, 37293, 37309, 37325, 37341, 37357, 37377, 37393, 37409, 37513, 37529, 37545, 37561, 37577, 37597, 37613, 37629, 37689, 37705, 37721, 37737, 37753, 37773, 37789, 37849]
                xData = [887.59, 889.13, 889.71, 890.29, 890.92, 891.59, 892.27, 894.78, 895.5, 896.22, 896.95, 897.69, 898.51, 899.29, 899.92, 902.1, 902.74, 903.37, 903.98, 904.58, 905.16, 905.7, 906.19, 907.75, 908.14, 908.47, 908.75, 909.03, 909.29, 909.5, 910.21, 910.23, 910.19, 910.13, 910.01, 909.81, 909.59, 909.38, 908.27, 907.95, 907.63, 907.23, 906.74, 906.24, 905.79, 904.17, 903.63, 903.04, 902.39, 901.71, 901.1, 900.47, 899.78, 897.25, 896.48, 895.71, 894.94, 894.18, 893.42, 892.63, 888.97, 888.33, 887.69, 887.06, 886.45, 885.88, 885.39, 884.91, 883.48, 883.16, 882.84, 882.57, 882.36, 882.23, 882.15, 882.08, 881.96, 881.99, 882.13, 882.31, 882.48, 882.68, 882.87, 883.82, 884.23, 884.68, 885.14, 885.63, 886.13, 886.65, 887.19, 890.34, 891.01, 891.75, 892.55, 893.35, 894.14, 894.94, 896.39, 898.01, 898.79, 899.56, 900.32, 901.05, 901.8, 902.54, 904.97, 905.53, 906.04, 906.5, 906.89, 907.31, 907.76, 908.13, 909.11, 909.33, 909.49, 909.6, 909.69, 909.74, 909.72, 909.25, 909.13, 908.98, 908.75, 908.48, 908.18, 907.89, 907.58, 906.27, 905.85, 905.39, 904.88, 904.39, 903.86, 903.32, 901.4, 900.8, 900.13, 899.43, 898.78, 898.14, 897.5, 896.82, 894.41, 893.76, 893.12, 892.48, 891.85, 891.17, 890.49, 889.84, 886.85, 886.33, 885.85, 885.4, 884.95, 884.48, 884.01, 882.74, 882.48, 882.23, 882.03, 881.87, 881.7, 881.53, 881.44, 881.36, 881.36, 881.42, 881.55, 881.72, 881.89, 882.06, 883.09, 883.45, 883.83, 884.43, 884.98, 885.72, 886.35, 886.94, 890.46, 891.23, 891.97, 892.68, 893.45, 894.29, 895.16, 896.04, 899.26, 900.29, 901.1, 901.88, 902.62, 903.3, 903.95, 906.16, 906.66, 907.11, 907.53, 907.91, 908.32, 908.62, 908.92, 909.63, 909.73, 909.74, 909.7, 909.68, 909.67, 909.62, 909.46, 908.19, 907.91, 907.64, 907.3, 906.94, 906.64, 906.29, 904.66, 904.21, 903.61, 903.11, 902.59, 902.04, 901.47, 900.89, 898.83, 898.25, 897.66, 896.93, 896.36, 895.78, 895.17, 892.83, 892.24, 891.69, 891.13, 890.52, 889.72, 889.11, 888.53, 885.11, 884.67, 884.28, 883.88, 883.49, 883.08, 882.74, 882.42, 881.65, 881.51, 881.38, 881.31, 881.27, 881.28, 881.32, 881.75, 881.16]
                yData = [1586.86, 1587.98, 1588.36, 1588.68, 1589, 1589.33, 1589.65, 1590.42, 1590.58, 1590.69, 1590.77, 1590.85, 1590.91, 1590.88, 1590.8, 1590.45, 1590.3, 1590.09, 1589.86, 1589.56, 1589.22, 1588.9, 1588.58, 1587.2, 1586.76, 1586.24, 1585.67, 1585.17, 1584.69, 1584.12, 1580.57, 1579.91, 1579.24, 1578.54, 1577.81, 1577.08, 1576.37, 1575.7, 1573.32, 1572.64, 1571.95, 1571.27, 1570.54, 1569.85, 1569.19, 1567.16, 1566.65, 1566.15, 1565.65, 1565.18, 1564.74, 1564.34, 1563.96, 1562.9, 1562.71, 1562.58, 1562.46, 1562.37, 1562.32, 1562.28, 1562.83, 1563.04, 1563.27, 1563.55, 1563.85, 1564.18, 1564.57, 1564.97, 1566.67, 1567.26, 1567.86, 1568.43, 1569.04, 1569.67, 1570.34, 1571.08, 1573.56, 1574.22, 1574.9, 1575.69, 1576.44, 1577.15, 1577.85, 1580.1, 1580.72, 1581.35, 1581.96, 1582.56, 1583.19, 1583.75, 1584.3, 1586.89, 1587.3, 1587.66, 1587.97, 1588.28, 1588.58, 1588.85, 1589.2, 1589.44, 1589.49, 1589.44, 1589.33, 1589.18, 1588.97, 1588.75, 1587.81, 1587.52, 1587.15, 1586.71, 1586.19, 1585.66, 1585.17, 1584.7, 1582.82, 1582.22, 1581.61, 1581.01, 1580.42, 1579.8, 1579.15, 1575.98, 1575.36, 1574.75, 1574.14, 1573.54, 1572.96, 1572.38, 1571.81, 1569.68, 1569.09, 1568.52, 1567.99, 1567.52, 1567.02, 1566.5, 1564.87, 1564.47, 1564.16, 1563.89, 1563.54, 1563.15, 1562.83, 1562.58, 1562.06, 1561.96, 1561.89, 1561.84, 1561.78, 1561.75, 1561.81, 1561.93, 1562.98, 1563.29, 1563.62, 1563.94, 1564.28, 1564.64, 1565.02, 1566.62, 1567.08, 1567.55, 1568.05, 1568.58, 1569.14, 1569.68, 1570.18, 1572.16, 1572.79, 1573.48, 1574.17, 1574.85, 1575.52, 1576.18, 1578.69, 1579.38, 1580.11, 1581.01, 1581.7, 1582.52, 1583.15, 1583.72, 1586.25, 1586.71, 1587.18, 1587.59, 1587.97, 1588.29, 1588.54, 1588.75, 1589.03, 1588.98, 1588.89, 1588.76, 1588.61, 1588.43, 1588.2, 1586.94, 1586.51, 1586.08, 1585.67, 1585.27, 1584.71, 1584.17, 1583.53, 1581.21, 1580.6, 1579.94, 1579.27, 1578.61, 1577.96, 1577.31, 1576.5, 1573.1, 1572.49, 1571.93, 1571.34, 1570.65, 1569.96, 1569.34, 1567.43, 1566.95, 1566.37, 1565.97, 1565.57, 1565.15, 1564.73, 1564.34, 1563.05, 1562.75, 1562.5, 1562.25, 1562.04, 1561.85, 1561.7, 1561.32, 1561.29, 1561.28, 1561.3, 1561.35, 1561.48, 1561.58, 1561.71, 1563.58, 1564, 1564.44, 1564.89, 1565.37, 1565.9, 1566.32, 1566.79, 1568.82, 1569.37, 1569.94, 1570.58, 1571.21, 1571.97, 1572.6, 1575.1, 1566.09]
                zData = [424.65, 424.91, 425.05, 425.15, 425.24, 425.37, 425.57, 426.47, 426.72, 426.93, 427.14, 427.37, 427.62, 427.9, 428.19, 429.33, 429.66, 429.95, 430.25, 430.57, 430.88, 431.19, 431.54, 432.76, 433.13, 433.54, 433.92, 434.26, 434.56, 434.86, 436.42, 436.7, 436.97, 437.25, 437.52, 437.77, 437.99, 438.19, 438.85, 439.02, 439.16, 439.28, 439.38, 439.46, 439.57, 439.72, 439.7, 439.67, 439.6, 439.47, 439.38, 439.35, 439.31, 438.83, 438.62, 438.39, 438.15, 437.87, 437.58, 437.27, 435.39, 435.02, 434.7, 434.36, 434, 433.62, 433.21, 432.78, 431.34, 430.96, 430.6, 430.24, 429.84, 429.48, 429.16, 428.84, 427.79, 427.51, 427.26, 427, 426.75, 426.57, 426.44, 426.01, 425.91, 425.82, 425.7, 425.61, 425.55, 425.49, 425.45, 425.69, 425.8, 425.95, 426.12, 426.33, 426.55, 426.74, 427.16, 427.75, 428.05, 428.39, 428.81, 429.2, 429.55, 429.91, 431.31, 431.67, 432.02, 432.38, 432.72, 433.09, 433.45, 433.81, 435.05, 435.36, 435.68, 436.02, 436.32, 436.61, 436.92, 438.33, 438.53, 438.77, 438.97, 439.11, 439.29, 439.48, 439.61, 439.84, 439.93, 439.97, 439.96, 439.96, 439.95, 439.97, 439.8, 439.67, 439.55, 439.47, 439.33, 439.2, 439.08, 438.92, 438.25, 438.04, 437.84, 437.57, 437.28, 437, 436.68, 436.35, 434.66, 434.3, 433.91, 433.56, 433.23, 432.91, 432.58, 431.48, 431.21, 430.9, 430.59, 430.26, 429.91, 429.61, 429.34, 428.46, 428.2, 427.98, 427.74, 427.5, 427.33, 427.16, 426.47, 426.32, 426.18, 426, 425.91, 425.8, 425.74, 425.75, 426.1, 426.22, 426.34, 426.47, 426.66, 426.88, 427.1, 427.35, 428.54, 428.98, 429.35, 429.7, 430.04, 430.41, 430.79, 432.22, 432.61, 432.99, 433.35, 433.71, 434.16, 434.5, 434.81, 436.04, 436.39, 436.7, 436.95, 437.26, 437.6, 437.9, 438.22, 439.05, 439.15, 439.24, 439.37, 439.52, 439.68, 439.78, 439.89, 439.95, 440.02, 440.02, 439.95, 439.88, 439.85, 439.85, 439.58, 439.43, 439.29, 439.12, 438.97, 438.8, 438.59, 437.8, 437.55, 437.27, 437.06, 436.82, 436.41, 436.1, 435.79, 433.47, 433.13, 432.81, 432.47, 432.13, 431.72, 431.38, 431.05, 429.9, 429.59, 429.32, 429.04, 428.72, 428.38, 428.13, 427.27, 430.74]
                xTest = [881.04, 880.44, 889.41, 901.18, 910.09, 911.06, 901.67, 889.44, 880.74, 880.39, 888.66, 900.99, 910.05, 910.09, 901.37, 889.26, 880.29, 879.88, 888.69, 900.51, 909.99, 909.77, 901.22, 888.58, 879.79]
                yTest = [1566.06, 1577.64, 1588.18, 1591.76, 1585.86, 1573.95, 1562.89, 1559.31, 1565.21, 1576.89, 1587.86, 1591.16, 1584.79, 1573.18, 1562.72, 1559.14, 1565.25, 1576.5, 1587.41, 1590.28, 1584.86, 1572.92, 1562.46, 1558.41, 1564.11]
                zTest = [430.54, 425.51, 424.33, 428.21, 434.21, 439.48, 440.53, 437.04, 430.59, 425.85, 424.91, 428.08, 434.46, 439.11, 440.31, 436.91, 430.22, 425.46, 424.44, 427.63, 434.36, 439.48, 440.25, 436.84, 430.5]
                break
        }

        // transpose the three arrays into array of triples
        for (let i = 0; i < scanTimes.length; i++) {
            let xyz = []
            xyz.push(xData[i])
            xyz.push(yData[i])
            xyz.push(zData[i])
            scanData.push(xyz)
            // when simulating, only bother Tracing every 10th sample
            if ((i % 10) == 0) {
                datalogger.log(
                    datalogger.createCV("t", scanTimes[i]),
                    datalogger.createCV("x", xData[i]),
                    datalogger.createCV("y", yData[i]),
                    datalogger.createCV("z", zData[i])
                )
            }
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
