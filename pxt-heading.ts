/**
 * An extension providing a compass-bearing for a buggy located anywhere on the globe 
 * (except at the magnetic poles!), with any arbitrary mounting orientation for its microbit.
 * See the README for a detailed description of the approach, methods and algorithms.
 * 
 * TODO? No use is yet made of the accelerometer. Although helpful to compensate for static tilt,
 * on a moving buggy dynamic sideways accelerations confound the measurement of "down", so 
 * applying tilt-compensation could actually make compass-heading accuracy worse!
 */

// OPERATIONAL MODES (controlling data-logging)
enum Mode {
    Normal, // Normal usage, mounted in a buggy
    Capture, // Acquire a new test dataset, using a special rotating jig
    Debug, // Test & debug (NOTE: test-dataset selection is hard-coded below)
}

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
    const ThreePi = 3 * Math.PI
    const HalfPi = Math.PI / 2
    const RadianDegrees = 360 / TwoPi
    const EnoughScanTime = 1800 // minimum acceptable scan-time
    const EnoughSamples = 70 // fewest acceptable scan samples
    const MarginalField = 30 // minimum acceptable field-strength for magnetometer readings
    const Circular = 1.03 // maximum eccentricity to consider an Ellipse as "circular"
    const LongEnough = 0.9 // for major-axis candidates, qualifying fraction of longest 

    // SUPPORTING CLASSES

    // An Arrow is an object holding a directed vector {u,v} in both Cartesian and Polar coordinates. 
    // It also carries a time field, used to timestamp scanned samples.
    // It is used to hold a 2D magnetometer measurement as a normalised vector.

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
    // view formed when projecting the scan Spin-Circle onto a 2-axis View-plane.
    class Ellipse {
        plane: string; // View name (just for debug)
        uDim: number; // horizontal axis of this View
        vDim: number; // vertical axis of this View
        uOff: number; // horizontal offset needed to re-centre this Ellipse along the U-axis
        vOff: number; // vertical offset needed to re-centre this Ellipse along the V-axise circular again

        // calibration characteristics
        majorAxis: Arrow; // averaged major axis 
        period: number; // this View's assessment of average scan-rotation time
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

        /* Method to find the average of a set of adjacent Arrow angles (coping with cyclic angle roll-round)
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
        */

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
                // the first candidate fixes which "end" of the axis we're choosing
                let front = arrows[0].angle
                uSum = arrows[0].u
                vSum = arrows[0].v
                let startTime = arrows[0].time
                for (let i = 1; i < count; i++) {
                    // get angle difference, as +/- 2pi
                    let deviate = ((ThreePi + arrows[i].angle - front) % TwoPi) - Math.PI
                    // does it point mostly to the front? ...or to the back?
                    if (Math.abs(deviate) < HalfPi) {
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
                    if (mode == Mode.Debug) {
                        datalogger.log(
                            datalogger.createCV("time", arrows[i].time),
                            datalogger.createCV("angle", round2(arrows[i].angle)),
                            datalogger.createCV("bearing", round2(asBearing(arrows[i].angle))),
                            datalogger.createCV("uSum", round2(uSum)),
                            datalogger.createCV("vSum", round2(vSum)),
                            datalogger.createCV("turns", turns))
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
        // It also collects possible candidates for the Ellipse major-axis by looking for inflections in the 
        // slope that occur as we pass local peaks; candidate values are pushed onto this.majors[] axis-list
        // for later averaging and timing analysis.
        extractAxes() {
            let majors: Arrow[] = [] // candidate directions for major axis of Ellipse
            let peak = 0 // (a marker, just for debug trace)

            let trial = new Arrow(scanData[0][this.uDim], scanData[0][this.vDim], scanTimes[0])
            let longest = trial.size
            let shortest = trial.size
            if (mode == Mode.Debug) {
                datalogger.log(
                    datalogger.createCV("time", trial.time),
                    datalogger.createCV("size", round2(trial.size)),
                    datalogger.createCV("bearing", round2(asBearing(trial.angle))),
                    datalogger.createCV("longest", round2(longest)))
            }
            let sizeWas: number
            let step: number = 99999 // marker for "first time round"
            let stepWas: number
            // compare first with remaining samples
            for (let i = 1; i < scanTimes.length; i++) {
                sizeWas = trial.size
                trial = new Arrow(scanData[i][this.uDim], scanData[i][this.vDim], scanTimes[i])
                longest = Math.max(longest, trial.size)
                shortest = Math.min(shortest, trial.size)
                // now collect candidates for the major-axis angle   
                stepWas = step
                step = trial.size - sizeWas // is radius growing or shrinking?
                // ensure the first two steps will always match
                if (stepWas == 99999) stepWas = step
                // look for peaks where we switch from growing to shrinking
                if ((stepWas > 0) && (step <= 0)) {
                    majors.push(trial.cloneMe()) // copy the major axis we are passing
                    if (mode == Mode.Debug) {
                        datalogger.log(
                            datalogger.createCV("time", trial.time),
                            datalogger.createCV("size", round2(trial.size)),
                            datalogger.createCV("bearing", round2(asBearing(trial.angle))),
                            datalogger.createCV("longest", round2(longest)))
                    }
                }
            }

            // The ratio of the axis lengths gives the eccentricity of this Ellipse
            this.eccentricity = longest / shortest
            // Readings taken from a near-circular Ellipse won't be fore-shortened, so we can skip correction!
            this.isCircular = (this.eccentricity < Circular)

            /* We are trying to find a good approximation to the tilt of the Ellipse's major-axis.  
            We could simply nominate the longest one detected, but instead we average the candidates. 
            Passing the major-axis twice per Spin-circle revolution, an eccentric Ellipse will produce neatly
            alternating candidates with "opposite" angles. Noisy readings mean that a more-nearly circular
            Ellipse may generate alternating clusters of candidates as we pass each end of the axis.
            An almost circular Ellipse has no meaningful axis, and will generate multiple spurious candidates. 
            */
            // purge any local maximum whose vector length is too short --it's nowhere near the major-axis
            let long = longest * LongEnough
            for (let i = 0; i < majors.length; i++) {
                if (majors[i].size < long) {
                    majors.splice(i, 1)  // disqualified!
                    i-- // (all subsequent candidates shuffle up by one place!)
                }
            }
            this.majorAxis = this.combine(majors)  // form a consensus 
            this.period = this.majorAxis.time 

            if (mode == Mode.Debug) {
                let theta = this.majorAxis.angle
                datalogger.log(
                    datalogger.createCV("view", this.plane),
                    datalogger.createCV("theta", round2(theta)),
                    datalogger.createCV("bearing", round2(asBearing(theta))),
                    datalogger.createCV("period", round2(this.period)),
                    datalogger.createCV("eccent.", round2(this.eccentricity)))
            }
        }
        // By comparing the longest and shortest radii, this method works out the eccentricity of this Ellipse.
        // It also collects possible candidates for the Ellipse major-axis into this.majors[]
        // XXX Although samples are already a rolling sum of seven readings, differentiation will amplify noise,
        // XXX so for stability we apply an additional 7-point Gaussian smoothing as we track the first-derivative.
        // We look for inflections in the slope, which occur as we pass local peaks; then push candidate values 
        // onto the majors[] axis-list for later averaging and timing analysis.
        extractAxesWas() {
            let arrows: Arrow[] = [] // rolling set of 7 Arrows, holding 7 adjacent samples
            let majors: Arrow[] = [] // candidate directions for major axis of Ellipse
            let peak = 0 // (a marker, just for debug trace)
            // collect 7 samples as Arrows and get first Gaussian sum
            for (let i = 0; i < 7; i++) {
                arrows.push(new Arrow(scanData[i][this.uDim], scanData[i][this.vDim], scanTimes[i]))
            }
            // let smooth = this.gaussianAverage(arrows)
            let smooth = arrows[3] // try unsmoothed

            let longest = smooth.size
            let shortest = smooth.size
            let sizeWas: number
            let step: number = 99999 // marker for "first time round"
            let stepWas: number
            // "smooth" contains a smoothed (w.r.t its neighbours) version of the 4th sample in this View
            // now work through remaining samples...
            for (let i = 7; i < scanTimes.length; i++) {
                sizeWas = smooth.size
                arrows.shift() // drop the earliest sample & append the next one
                arrows.push(new Arrow(scanData[i][this.uDim], scanData[i][this.vDim], scanTimes[i]))
                // smooth = this.gaussianAverage(arrows)
                smooth = arrows[3]
                longest = Math.max(longest, smooth.size)
                shortest = Math.min(shortest, smooth.size)
                // now collect candidates for the major-axis angle   
                stepWas = step
                // is radius growing or shrinking?
                step = smooth.size - sizeWas
                // ensure the first two steps always match
                if (stepWas == 99999) stepWas = step
                // look for peaks where we switch from growing to shrinking
                if ((stepWas > 0) && (step <= 0)) {
                    majors.push(smooth.cloneMe()) // copy the major axis we are passing
                    if (mode == Mode.Debug) {
                        datalogger.log(
                            datalogger.createCV("time", smooth.time),
                            datalogger.createCV("size", round2(smooth.size)),
                            datalogger.createCV("bearing", round2(asBearing(smooth.angle))),
                            datalogger.createCV("longest", round2(longest)))
                    }
                }
            }
            // The ratio of the axis lengths gives the eccentricity of this Ellipse
            this.eccentricity = longest / shortest
            // Readings taken from a near-circular Ellipse won't be fore-shortened, so we can skip correction!
            this.isCircular = (this.eccentricity < Circular)

            /* We are trying to find a good approximation to the tilt of the Ellipse's major-axis.  
            We could simply nominate the longest one detected, but instead we average the candidates. 
            Passing the major-axis twice per Spin-circle revolution, an eccentric Ellipse will produce neatly
            alternating candidates with "opposite" angles. Noisy readings mean that a more-nearly circular
            Ellipse may generate alternating clusters of candidates as we pass each end of the axis.
            An almost circular Ellipse has no meaningful axes, and will generate multiple spurious candidates. 
            */
            // purge any local maximum whose vector length is too short --it's nowhere near the major-axis
            let long = longest * LongEnough
            for (let i = 0; i < majors.length; i++) {
                if (majors[i].size < long) {
                    majors.splice(i, 1)  // disqualified!
                    i-- // (all subsequent candidates shuffle up by one place!)
                }
            }

            this.majorAxis = this.combine(majors)
            this.period = this.majorAxis.time


            if (mode == Mode.Debug) {
                let theta = this.majorAxis.angle
                datalogger.log(
                    datalogger.createCV("view", this.plane),
                    datalogger.createCV("theta", round2(theta)),
                    datalogger.createCV("bearing", round2(asBearing(theta))),
                    datalogger.createCV("period", round2(this.period)),
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

    let mode:Mode = Mode.Debug // mode switch for logging
    //let logging = true  // logging mode flag
    //let capturing = false // data capturing mode
    //let debugging = false  // test mode flag
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

        if (mode != Mode.Normal) {
            datalogger.deleteLog()
            datalogger.includeTimestamp(FlashLogTimeStampFormat.None)
        }

        if (mode == Mode.Debug) {
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

        if (mode == Mode.Capture) {
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
        
/*        bestView = View.ZX  // while debugging, force use of this view!
*/
        // periodicity is unreliable in the best View: average just the other two Views' measurements
        period = (views[0].period + views[1].period + views[2].period - views[bestView].period) / 2
        rpm = 60000 / period

        // Depending on mounting orientation, the bestView Ellipse might possibly be seeing the 
        // Spin-Circle from "underneath", effectively experiencing an anti-clockwise scan.
        // Check polarity of revolution:
        fromBelow = (views[bestView].period < -1)

        // For efficiency, extract various characteristics from the bestView Ellipse
        uDim = views[bestView].uDim
        vDim = views[bestView].vDim
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

        if (mode == Mode.Debug) {
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
     * Choose mode: whether to run normally,
     *  - or to use Data Logger to grab a new test dataset into MY_DATA,
     *  - or to debug processing using a named test dataset
     */
    //% block="set mode $mode"
    //% inlineInputMode=inline 
    //% weight=10
    export function setMode(newMode:Mode, name: string) {
        mode = newMode
        dataset = name
    }




    // UTILITY FUNCTIONS
    // Take the sum of seven new readings to get a stable fix on the current heading

    /* Although eventually we'd only need [uDim, vDim], we sum and log all three dims.
       This will allow us, while testing, to override automatic choice of bestView
       and check out more severe levels of correction!
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

            if (mode == Mode.Capture) {
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
            // Now scale up along V, re-balancing the axes to make the Ellipse circular
            vFix = vNew * scale
            // get the adjusted angle for this corrected {u,v}
            reading = Math.atan2(vFix, uFix)
            // finally, undo the rotation by theta
            reading += theta
        } else {
            reading = Math.atan2(v, u)
        }

        if (mode == Mode.Debug) {
            datalogger.log(
                datalogger.createCV("u", round2(u)),
                datalogger.createCV("v", round2(v)),
                datalogger.createCV("uNew", round2(uNew)),
                datalogger.createCV("vNew", round2(vNew)),
                datalogger.createCV("uFix", round2(uFix)),
                datalogger.createCV("vFix", round2(vFix)),
                datalogger.createCV("reading", round2(reading)),
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
