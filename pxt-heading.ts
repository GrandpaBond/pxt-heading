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
    const Window = 7 // number of magnetometer samples needed to form a good average
    const SampleGap = 15 // minimum ms to wait between magnetometer readings
    const AverageGap = 25 // (achieved in practice, due to system interrupts)
    const Latency = Window * AverageGap // resulting time taken to collect a good moving average from scratch

    // GLOBALS

    // (For testing purposes scan data is made externally visible)
    export let scanTimes: number[] = [] // sequence of time-stamps for scanned readings 
    export let scanData: number[][] = [] // scanned sequence of [X,Y,Z] magnetometer readings
    export let nSamples: number

    // (Other externally visible test-related globals)
    export let debugMode = false // in debugMode we use pre-loaded data
    export let test = 0 // global selector for single readings
    export let testTimes: number[] = [] // sequence of time-stamps for single readings
    export let testData: number[][] = [] //[X,Y,Z] magnetometer values for single readings

    // (while debugging, make the Ellipse data externally visible too)
    export let xy: Ellipse
    export let yz: Ellipse
    export let zx: Ellipse


    let plane = "undefined"
    let uDim = -1 // the "horizontal" axis (called U) for the best View
    let vDim = -1 // the "vertical" axis (called V) for the best View
    let north = 0 // reading registered as "North"
    let strength = 0 // the average magnetic field-strength observed by the magnetometer
    let period = -1 // overall assessment of average rotation time

    // amplitudes and central offsets of sinusoidal scan-readings in each dimension
    let xField: number
    let yField: number
    let zField: number
    let xOff: number
    let yOff: number
    let zOff: number

    // correction parameters adopted from bestView Ellipse for future readings
    let rotationSense = 1 // set to -1 if orientation means field-vector projection is "from below"
    let isCircular: boolean // if best view Ellipse is circular, no correction is needed
    let uOff: number // horizontal origin offset
    let vOff: number // vertical origin offset
    let theta: number  // major-axis tilt angle (in radians anticlockwise from the U-axis)
    let cosTheta: number; // saved for efficiency
    let sinTheta: number; //      ditto
    let scale: number // stretch-factor for correcting foreshortened readings (= eccentricity)


    // EXPORTED USER INTERFACES   

    /** 
     * Assuming the buggy is currently spinning clockwise on the spot, capture a time-stamped
     * sequence of magnetometer readings, and when finished, process them to set up the compass.
     *
     * @param ms scanning-time in millisecs (long enough for more than one full rotation)
     * @return zero if successful, or a negative error code:
     *
     *      -1 : NOT ENOUGH SCAN DATA

     *      -2 : FIELD STRENGTH TOO WEAK
     *
     *      -3 : NOT ENOUGH SCAN ROTATION
     *
     */

    //% block="scan clockwise for (ms) $ms" 
    //% inlineInputMode=inline 
    //% ms.shadow="timePicker" 
    //% ms.defl=0 
    //% weight=90 

    export function scanClockwise(ms: number): number {
        strength = -1
        period = -1

        // unless test-data has already been pre-loaded...
        if (!debugMode) collectSamples(ms)  // ...take repeated magnetometer readings
        
        nSamples = Math.min(scanTimes.length, scanData.length) // (defend against pre-load mismatch)

        // Now analyse the scan-data to decide how best to use the magnetometer readings.
        // we'll typically need about a couple of second's worth of scanned readings...
        let scanTime = scanTimes[nSamples - 1] - scanTimes[0]
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
        xField = (xhi - xlo) / 2
        yField = (yhi - ylo) / 2
        zField = (zhi - zlo) / 2
        strength = Math.sqrt((xField * xField) + (yField * yField) + (zField * zField))

        // Bail out early if the scan didn't properly detect the Earth's magnetic field,
        // (perhaps due to magnetic shielding?)
        if (strength < MarginalField) {
            return -2 // "FIELD STRENGTH TOO WEAK"
        }

        // The means of the extremes give a good approximation to the central offsets.
        xOff = (xhi + xlo) / 2
        yOff = (yhi + ylo) / 2
        zOff = (zhi + zlo) / 2

        // 2nd pass re-centres all the scanData samples, eliminating "hard-iron" environmental magnetic effects.
        for (let i = 0; i < nSamples; i++) {
            scanData[i][Dimension.X] -= xOff
            scanData[i][Dimension.Y] -= yOff
            scanData[i][Dimension.Z] -= zOff
        }

        // create three Ellipse instances, for analysing each possible 2D view of the spin-Circle
        xy = new Ellipse("XY", Dimension.X, Dimension.Y, xOff, yOff)
        yz = new Ellipse("YZ", Dimension.Y, Dimension.Z, yOff, zOff)
        zx = new Ellipse("ZX", Dimension.Z, Dimension.X, zOff, xOff)

        extractAxes() // do as it says on the tin...

        // check that at least one view saw at least one complete rotation (with a measurable period)...
        if ((xy.period + yz.period + zx.period) < 0) {
            return -3 // "NOT ENOUGH SCAN ROTATION"
        }

        // Choose the "roundest" Ellipse view  --the one with lowest eccentricity.
        let view: Ellipse = xy
        if (yz.eccentricity < view.eccentricity) view = yz
        if (zx.eccentricity < view.eccentricity) view = zx

        // periodicity may be unreliable in a near-circular View, so we form an average using
        // just the other two views' measurements
        period = (xy.period + yz.period + zx.period - view.period) / 2

        // For efficiency, extract various characteristics from the best Ellipse
        plane = view.plane
        uDim = view.uDim
        vDim = view.vDim
        uOff = view.uOff
        vOff = view.vOff
        rotationSense = view.rotationSense // whether reading from above or below the plane
        scale = view.eccentricity // the scaling needed to balance axes
        isCircular = (scale <= Circular)
        theta = view.tilt // the rotation (in radians) of the major-axis from the U-axis
        cosTheta = view.cosa
        sinTheta = view.sina

        // NOTE: although we have now finished with the scanData, its memory is not released yet
        //       but left accessible for potential capture as a test dataset.

        return 0
    }

    /**
     * Read the magnetometer and register the buggy's current direction as "North",
     * (i.e. the direction that will in future return zero as its heading).
     * 
     * The actual direction of the buggy when this function is called is arbitrary:
     * it could be Magnetic North; or True North (compensating for local declination); 
     * or any convenient direction from which to measure subsequent heading angles.
     * 
     * @return zero if successful, or a negative error code:
     *
     *      -4 // ERROR: SUCCESSFUL SCAN IS NEEDED FIRST
     */
    //% block="set North" 
    //% inlineInputMode=inline 
    //% weight=80 
    export function setNorth(): number {
        test = 0 // reset test history

        if (period == -1) {
            return -4 // ERROR: SUCCESSFUL SCAN IS NEEDED FIRST
        } else {
            // Having successfully set up the projection parameters for the bestView, get a
            // stable fix on the current heading, which we will then designate as "North".
            // (This is the global fixed bias to be subtracted from all future readings)
            north = 0
            north = takeSingleReading()
        }

        if (!debugMode) {
            // we've now definitely finished with the scanning data, so release its memory
            scanTimes = []
            scanData = []
        }
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

        if (period == -1) {
            return -4 // ERROR: SUCCESSFUL SCAN IS NEEDED FIRST
        } else {
            let heading = (takeSingleReading() - north) * rotationSense

            /* NOTE: that there is a double reversal going on here:
            Viewed from above, the Field-vector reading in radians INCREASES (anticlockwise) w.r.t "North"
            as the buggy's compass-heading INCREASES (clockwise).
            The reading will therefore increase by HalfPi after a right-turn from initially facing "North".
            So subtracting North (cyclically), that converts asDegrees() to +90.
            From below (when rotationSense = -1), the same right-turn would DECREASE the reading by HalfPi
            necessitating a third reversal, but only after first subtracting North !
            */

            return asDegrees(heading)
        }
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

    // SUPPORTING CLASSES

    // A Smoother object computes a moving average from a sequence of time-stamped values: 
    // in our case, magnetometer readings and their derivatives.
    // Timing irregularites due to scheduler interrupts demand this somewhat complex maths.
    // The constant {Window} governs the latency of the exponential averaging process.
    // Smoothers can work with arbitrary-sized vectors of values that share the same timestamp.
    // history[], previous[], latest[] and result[] arrays will be either 3-D, (for initial [X,Y,Z] scanning)
    // or 2-D, (for analysing the chosen [uDim,vDim]).
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
    An Ellipse is an object holding the characteristics of the view formed when projecting the 
    magnetic field-vector onto a 2-axis View-plane {XY, YZ or ZX}.
    
    While scanning clockwise, the projection of the field-vector will appear to trace out a (typically)
    elliptical view of the Spin-Circle, anti-clockwise. If viewed from "above" the polar angles of
    successive readings (in radians) will INCREASE (anti-clockwise); if viewed from "below" they DECREASE.

    The foreshortened view means that evenly spaced heading angles will appear bunched around the ends 
    of the ellipse. (In the extreme "side-on" case, the ellipse collapses to just a straight line!)

    To correct for foreshortening we must stretch points on the ellipse back onto a circle, by rescaling
    parallel to the minor axis of the ellipse until it matches the major axis. 
    The ratio of major-axis to minor-axis gives the necessary scaling factor {eccentricity}.

    Depending on the exact mounting orientation of the microbit in the buggy, this ellipse may be
    tilted with respect to a particular view's axes, so correction then becomes a three stage process:
    1) rotate the new point [u,v] by {-tilt} so the minor-axis lines up with the V-axis.
    2) stretch the V coordinate by {eccentricity}
    3) rotate back by {tilt} to give the corrected [u,v] from which the heading can be derived.
  
    This correction can in theory be applied in any of the three views (unless exactly side on to
    the Spin-Circle), but the best accuracy is obtained from the most circular of the three views.
    Readings on a near-circular Ellipse are barely fore-shortened at all, so we can skip correction!
    */
    export class Ellipse {  // (make class definition externally visible while still debugging)
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

    
        newTurn: boolean  // allow clocking of a major-axis crossing
        turns: number // number of full rotations since the first major-axis transit
        start: number // timestamp of first major-axis transit
        finish: number // timestamp of latest major-axis transit

        // calibration characteristics
        period: number // scan-rotation time as viewed in this plane
        tilt: number // direction of the Ellipse major-axis 
        cosa: number // helpful cosine of the tilt
        sina: number // helpful sine of the tilt
        eccentricity: number // ratio of major-axis to minor-axis magnitudes for this Ellipse
        rotationSense: number // rotation sign = +/-1, reflecting this Ellipse's view of the clockwise scan


        constructor(plane: string, uDim: number, vDim: number, uOff: number, vOff: number) {
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
            this.tilt = 0
            this.start = 0
            this.finish = 0
            this.turns = -1
            this.newTurn = true
            this.period = -1
            this.rotationSense = 0
        }

        // addMajor() is called whenever we pass a major-axis (dQ changes from positive to negative)
        // It builds a resultant vector (reversing coordinates at the "other" end of the axis).
        addMajor(i: number, u: number, v: number, dw: number) {
            if (dw > 0) { // just surfacing above our plane so add-in current vector coordinates
                this.uHi += u
                this.vHi += v
                if (this.newTurn) { // clock this major-axis crossing
                    if (this.start == 0) {
                        this.start = scanTimes[i] // start measuring turns
                    }
                    this.finish = scanTimes[i]
                    this.turns++
                    this.newTurn = false // prevent future clocking of another major-axis crossing...
                    //... until after the next minor axis crossing
                }
            } else {  // just dipping below our plane so subtract current vector coordinates
                // (as it's the other end of the major-axis)
                this.uHi -= u
                this.vHi -= v
            }
            this.nHi++
        }

        // addMinor() is called whenever the Normal coordinate changes from growing to shrinking.
        // This means the field-vector is farthest from our plane (above or below), so we're near
        // the minor-axis of this plane's ellipse.
        addMinor(i: number, u: number, v: number, w: number) {
            if (w > 0) { // orthogonal reading is near its positive maximum
                this.uLo += u
                this.vLo += v
            } else { // orthogonal reading is near its negative maximum
                this.uLo -= u
                this.vLo -= v
            }
            this.nLo++
            this.newTurn = true // permit future clocking of a major-axis crossing
        }


        // calculate() method is called after all scandata has been processed...
        calculate() {
            this.eccentricity = -1
            this.period = -1
            this.tilt = 0
            if (this.nHi > 0) {
                // get average major-axis radius
                this.uHi /= this.nHi
                this.vHi /= this.nHi
                this.rHi = Math.sqrt(this.uHi * this.uHi + this.vHi * this.vHi)

                if (this.nLo > 0) {
                    // get average minor-axis radius
                    this.uLo /= this.nLo
                    this.vLo /= this.nLo
                    this.rLo = Math.sqrt(this.uLo * this.uLo + this.vLo * this.vLo)
                    // ratio of axes gives eccentricity
                    this.eccentricity = this.rHi / this.rLo

                    // use sign of  cross-product of axes to determine rotation-rotationSense
                    let cross = (this.uHi * this.vLo) - (this.uLo * this.vHi)
                    this.rotationSense = Math.abs(cross) / cross

                    // save the major-axis angle, and its components
                    this.tilt = Math.atan2(this.vHi, this.uHi)
                    this.cosa = this.uHi / this.rHi
                    this.sina = this.vHi / this.rHi
                }
            }
            // get average time for scan rotation 
            if (this.turns > 0) {
                this.period = (this.finish - this.start) / this.turns
            }

            /* logging of view characteristics used while testing only
            if(debugMode) {
                datalogger.log(
                    datalogger.createCV("plane", this.plane),
                    datalogger.createCV("uHi", round2(this.uHi)),
                    datalogger.createCV("vHi", round2(this.vHi)),
                    datalogger.createCV("uLo", round2(this.uLo)),
                    datalogger.createCV("vLo", round2(this.vLo)),
                    datalogger.createCV("tilt", round2(this.tilt)),
                    datalogger.createCV("eccen.", round2(this.eccentricity)),
                    datalogger.createCV("sense", this.rotationSense),
                    datalogger.createCV("turns", this.turns),
                    datalogger.createCV("period", round2(this.period))
                )
            }
            */
        }
    }


    // UTILITY FUNCTIONS

    /* Take magnetometer readings periodically over the specified duration (generally a couple
      of seconds), and append a sequence of new [X,Y,Z] triples to the scanData[] array.
      A timestamp for each sample is also recorded in the scanTimes[] array.

      NOTE: To smooth out jitter, each reading is always a moving average of several consecutive readings.
      Because sample-times may be irregular (due to scheduled interrupts), a Smoother is used to provide
      a timing-aware exponential moving-average. The sample-grouping and spacing are controlled 
      respectively by the constants Window and SampleGap, which together determine the Latency.

     NOTE:
      If we are in debug-mode the scanData[] and scanTimes[] will have been pre-loaded.

      This function is exported to allow new test datasets to be captured 
    */
    export function collectSamples( ms: number) {
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
            }
        }
    }

    /** Take the average of several new readings to get a stable fix on the current heading.
     *  @return the angle of the magnetic field-vector (in radians anticlockwise
     * from the horizontal U-axis), corrected for any fore-shortening due to projection
     * onto the bestView plane.
     * 
     * NOTE: leaves the latest reading in xyz[] to allow possible test-data capture
     */

    export function takeSingleReading(): number {
        let xyz: number[] = [0, 0, 0]
        let u = 0
        let v = 0
        let uNew = 0
        let vNew = 0
        let uFix = 0
        let vFix = 0
        let reading = 0

        if (debugMode) { // just choose the next test-data value (cyclically)
            xyz[Dimension.X] = testData[test][Dimension.X]
            xyz[Dimension.Y] = testData[test][Dimension.Y]
            xyz[Dimension.Z] = testData[test][Dimension.Z]
            test = (test + 1) % testData.length
        } else {
            // build a new sample as the average of {Window} consecutive 2D readings, {SampleGap} apart
            xyz = [0, 0, 0]
            for (let i = 0; i < Window; i++) {
                basic.pause(SampleGap)
                xyz[Dimension.X] += input.magneticForce(Dimension.X)
                xyz[Dimension.Y] += input.magneticForce(Dimension.Y)
                xyz[Dimension.Z] += input.magneticForce(Dimension.Z)
            }
            xyz[Dimension.X] /= Window
            xyz[Dimension.Y] /= Window
            xyz[Dimension.Z] /= Window

            // keep global history of single test readings (for possible later capture)
            testTimes.push(input.runningTime())
            testData.push(xyz)
            test++
        }

        // now pick the coordinates we actually want for the current view,
        // and re-centre this latest point w.r.t our Ellipse origin
        u = xyz[uDim] - uOff
        v = xyz[vDim] - vOff

        // Unless this Ellipse.isCircular, any {u,v} reading will be foreshortened in this View, and
        // must be stretched along the Ellipse minor-axis to place it correctly onto the Spin-Circle.

        if (isCircular) {
            reading = Math.atan2(v, u)
        } else {
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
            reading = (reading + theta + TwoPi) % TwoPi
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

    function extractAxes() {
        /* Each of our three views {XY, YZ, ZX} sees the Spin-circle as an ellipse.
        For correction of future readings we will need to find the eccentricity and axis-tilt of
        each ellipse, and then (for highest accuracy) choose the most circular view. 

        This function finds the major-axes (where the radius is at a maximum), and the minor-axes (where 
        the radius is at a minimum). To avoid multiple Math.sqrt() calls, we use the radius-squared 
        (the sum of the squares of the coordinates) as a proxy, called the Q-value.

        To find the axes we track the slope of Q, looking for inflections at peaks and troughs.
        Because of noisy data, fluctuations (or "bounces") can occur near an axis, especially for 
        more circular Ellipses. So a peak might occur near a minor-axis or a trough near a major-axis.
        At a maxor-axis the field is most nearly aligned with that plane, so the orthogonal reading 
        should always be near its smallest there. Conversely, at a minor-axis, the orthogonal reading should 
        always be near its peak amplitude. We use this fact to weed out spurious peaks and troughs.
        
        For better accuracy, multiple axis-crossings are averaged. Each axis gets passed TWICE per rotation,
        so it is important to either add or subtract vectors, according to which "end" is being crossed. 
        This is policed using the orthogonal reading (its slope for a major-axis; its sign for a minor-axis). 

        This function also calculates the apparent scan-rotation period by monitoring the times and count of 
        alternate major-axis crossings (ignoring any multiple contributions due to "bounces").

        The cross-product of its final axis-angles also gives us the rotation-sense as seen by each view.

        The function uses the global scanTimes[] and scanData[] arrays, and updates the three
        global Ellipse objects, {xy, yz and zx}.
        */

        // preload first samples
        let t = scanTimes[0]
        let x = scanData[0][Dimension.X]
        let y = scanData[0][Dimension.Y]
        let z = scanData[0][Dimension.Z]
        let xsq = x * x
        let ysq = y * y
        let zsq = z * z
        let qXY = xsq + ysq
        let qYZ = ysq + zsq
        let qZX = zsq + xsq

        // monitor the slope of the Q-curves
        let dqXY: number
        let dqYZ: number
        let dqZX: number

        // maintain some history
        let xWas: number
        let yWas: number
        let zWas: number
        let qXYWas: number
        let qYZWas: number
        let qZXWas: number
        let dqXYWas: number
        let dqYZWas: number
        let dqZXWas: number

        for (let i = 1; i < nSamples; i++) {
            // update history
            xWas = x
            yWas = y
            zWas = z
            qXYWas = qXY
            qYZWas = qYZ
            qZXWas = qZX
            dqXYWas = dqXY
            dqYZWas = dqYZ
            dqZXWas = dqZX

            // get next sample
            t = scanTimes[i]
            x = scanData[i][Dimension.X]
            y = scanData[i][Dimension.Y]
            z = scanData[i][Dimension.Z]
            xsq = x * x
            ysq = y * y
            zsq = z * z
            qXY = xsq + ysq
            qYZ = ysq + zsq
            qZX = zsq + xsq
            dqXY = qXY - qXYWas
            dqYZ = qYZ - qYZWas
            dqZX = qZX - qZXWas

            // to aid detection of sign-change in the slope, we doctor any values that are exactly zero
            if (dqXY == 0) dqXY = dqXYWas
            if (dqYZ == 0) dqYZ = dqYZWas
            if (dqZX == 0) dqZX = dqZXWas

            // Look for peaks and troughs in the Q-values for each plane.

            // A peak only qualifes as a major-axis when the orthogonal field is small: 
            // --the sign of its slope says which end is which
            if ((dqXY < 0) && (dqXYWas > 0)) {
                if (Math.abs(z) < zField * 0.25) {
                    xy.addMajor(i, x, y, z - zWas)
                }
            }

            if ((dqYZ < 0) && (dqYZWas > 0)) {
                if (Math.abs(x) < xField / 2) {
                    yz.addMajor(i, y, z, x - xWas)
                }
            }

            if ((dqZX < 0) && (dqZXWas > 0)) {
                if (Math.abs(y) < yField / 2) {
                    zx.addMajor(i, z, x, y - yWas)
                }
            }

            // A trough only qualifies as a minor-axis when orthogonal field is big: 
            // --its sign says which end is which
            if ((dqXY > 0) && (dqXYWas < 0)) {
                if (Math.abs(z) > zField / 2) {
                    xy.addMinor(i, x, y, z) 
                }
            }

            if ((dqYZ > 0) && (dqYZWas < 0)) {
                if (Math.abs(x) > xField / 2) {
                    yz.addMinor(i, y, z, x)
                }
            }

            if ((dqZX > 0) && (dqZXWas < 0)) {
                if (Math.abs(y) > yField / 2) {
                    zx.addMinor(i, z, x, y)
                }
            }


            /*** logging of major axes used while testing only
            if (debugMode) {
                datalogger.log(
                    datalogger.createCV("i", i),
                    datalogger.createCV("xy.u", round2(xy.uHi)),
                    datalogger.createCV("xy.v", round2(xy.vHi)),
                    datalogger.createCV("xy.n", round2(xy.nHi)),
                    datalogger.createCV("yz.u", round2(yz.uHi)),
                    datalogger.createCV("yz.v", round2(yz.vHi)),
                    datalogger.createCV("yz.n", round2(yz.nHi)),
                    datalogger.createCV("zx.u", round2(zx.uHi)),
                    datalogger.createCV("zx.v", round2(zx.vHi)),
                    datalogger.createCV("zx.n", round2(zx.nHi))
                )
            }
             ***/
        }
        /*** logging of minor-axes used while testing only
        if (debugMode) {
            datalogger.log(
                datalogger.createCV("i", -1),
                datalogger.createCV("xy.u", round2(xy.uLo)),
                datalogger.createCV("xy.v", round2(xy.vLo)),
                datalogger.createCV("xy.n", round2(xy.nLo)),
                datalogger.createCV("yz.u", round2(yz.uLo)),
                datalogger.createCV("yz.v", round2(yz.vLo)),
                datalogger.createCV("yz.n", round2(yz.nLo)),
                datalogger.createCV("zx.u", round2(zx.uLo)),
                datalogger.createCV("zx.v", round2(zx.vLo)),
                datalogger.createCV("zx.n", round2(zx.nLo))
            )
        }
        ***/

        // use the collected vector-sums to compute average axes, and thence the ellipse characteristics
        xy.calculate()
        yz.calculate()
        zx.calculate()
    }
}
