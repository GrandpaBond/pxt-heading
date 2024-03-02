/**
 * An extension providing a compass-heading for a buggy anywhere on the globe 
 * (except at the magnetic poles), whatever the mounting orientation of the microbit.
 */

//% color=#6080e0 weight=40 icon="\uf14e" block="Heading" 
namespace heading {
    enum Dim { // ...for brevity
        X = Dimension.X,
        Y = Dimension.Y,
        Z = Dimension.Z
    }

    // OBJECT CLASSES 
    // a Limit is a local minimum or maximum in a scanned array of readings 
    class Limit {
        time: number
        value: number
        constructor() {
            this.time = 0   // time-stamp when this Limit was hit
            this.value = 0  // signed value of local maximum amplitude
        }
    }


    // an Axis holds characteristics of one magnetometer axis
    class Axis {
        dim: Dim // the dimension this axis is describing (X=0;Y=1;Z=2)
        bias: number // the fixed offset to re-centre the scanned wave-form
        amp: number  // the average amplitude of the scanned wave-form
        limits: Limit[] // the array of local extremes detected
        nLimits: number // the count of local extremes detected
        time0: number   // the timestamp of the first limit in this axis
        limit0: number  // the normalised value of the first limit in this axis
        period: number  // the RPM of the scanned wave-form in this axis

        constructor(dim: Dim) {
            // dim should be redundant, as it always matches index as axes[] array is built
            this.dim = dim
            this.limits = []
        }

        // Method to characterise the scanData for each axis
        characteriseWas(timeStamp: number[], scanData: number[]) {
            // First we extract local maxima and minima from a sequence of scanned wave-form readings
            // For element [d] in [...a,b,c,d,e...], we multiply the averaged slope behind it
            // (given by [d]-[a]), with the slope ahead (given by [e]-[b]). 
            // Whenever the sign of the slopechanges (i.e this product is negative) we record a Limit: 
            // a duple object comprising the peak value [d], together with its time-stamp.
            this.nLimits = 0
            let then = timeStamp[0] - 201
            // loop through from [d] onwards
            for (let i = 3; i <= timeStamp.length - 3; i++) {
                let ahead = (scanData[i + 2] - scanData[i - 1])
                let behind = (scanData[i + 1] - scanData[i - 2])
                // An inflection-point is where the product will be negative.
                // (Ignore any too-close crossings, arising from excessive jitter)
                if (((ahead * behind) <= 0)
                    && ((timeStamp[i] - then) > 200)) {
                    let newLimit = new Limit
                    newLimit.value = scanData[i]
                    newLimit.time = timeStamp[i]
                    datalogger.log(datalogger.createCV("dim", this.dim),
                                    datalogger.createCV("time", newLimit.time),
                                    datalogger.createCV("value", newLimit.value))

                    then = timeStamp[i]
                    this.limits.push(newLimit)
                    this.nLimits++
                }
            }

            // Use the averages of paired limits to set the magnitude and offset.
            // Needs at least three of each (resulting from a complete spin)
            this.amp = 0
            this.bias = 0
            if (this.nLimits > 3) {
                // 1st pass to get offset 
                let sum = 0
                let n = 0
                // use balanced pairs (skipping first limit if length is odd) 
                for (let i = this.nLimits % 2; i < this.nLimits; i++) {
                    sum += this.limits[i].value
                    n++
                }
                this.bias = sum / n

                // 2nd pass uses offset to give amplitude  
                sum = 0
                for (let i = 0; i < this.nLimits; i++) {
                    sum += Math.abs(this.limits[i].value - this.bias)
                }
                this.amp = sum / this.nLimits
            }
            // note the time & value(normlised) of the first limit
            this.time0 = this.limits[0].time
            this.limit0 = this.limits[0].value - this.bias

            // work out the periodicity
            let spans = this.nLimits - 1
            this.period = 2 * (this.limits[spans].time - this.limits[0].time) / spans
        }
    }

    // GLOBALS
    const MarginalField = 100 // minimum acceptable amplitude for sum of magnetometer readings

    let scanTimes: number[] = [] // sequence of time-stamps for scanned readings 
    let scanData: number[][]= [] // scanned sequence of magnetometer readings  

    let axes: Axis[] = []
    axes.push(new Axis(Dim.X))
    axes.push(new Axis(Dim.Y))
    axes.push(new Axis(Dim.Z))

    let uDim = -1 // will hold the horizontal axis for projecting readings
    let vDim = -1 // will hold the vertical axis for projecting readings
    let uOff = 0 // the offset needed to re-centre the major (uDim) axis
    let vOff = 0 // the offset needed to re-centre the minor (vDim) axis
    let theta = 0 // the angle to rotate readings so the major (U) axis is horizontal
    let scale = 0 // the scaling to apply to minor (V) axis readings to balance axes
    let toNorth = 0 // the anglular bias to be added so that North = 0
    let uFlip = 1 // set to -1 if uDim polarity is inverted
    let vFlip = 1 // set to -1 if vDim polarity is inverted
    let plane = "??" // the projection plane we are using: "XY","YZ" or "ZX"
    let testing = false  // test mode flag
    let test = 0         // selector for test sample

    function dumpData() {
        datalogger.deleteLog()
        datalogger.includeTimestamp(FlashLogTimeStampFormat.None)
        datalogger.setColumnTitles("t", "x", "y", "z")
        for (let i = 0; i < scanTimes.length; i++) {
            datalogger.log(datalogger.createCV("t", scanTimes[i]),
                datalogger.createCV("x", scanData[Dim.X][i]),
                datalogger.createCV("y", scanData[Dim.Y][i]),
                datalogger.createCV("z", scanData[Dim.Z][i]))
        }
    }


    // USER INTERFACE

    // use some sample data, while debugging...
    export function simulateScan() {
        scanTimes = [5637, 5665, 5693, 5721, 5749, 5777, 5805, 5833, 5861, 5889, 5917, 5945, 5973, 6001, 6029, 6057, 6087, 6113, 6141, 6169, 6197, 6225, 6253, 6281, 6309, 6337, 6365, 6393, 6421, 6449, 6477, 6505, 6533, 6561, 6589, 6617, 6645, 6673, 6701, 6729, 6757, 6785, 6813, 6841, 6869, 6897, 6925, 6956, 6981, 7009, 7037, 7065, 7093, 7121, 7149, 7177, 7205, 7233, 7261, 7289, 7317, 7345, 7373, 7401, 7429, 7460, 7485, 7513, 7541, 7569, 7597, 7625, 7653, 7681, 7709, 7737, 7765, 7793, 7821, 7849, 7877, 7905, 7933, 7961, 7989, 8017, 8045, 8073, 8101, 8129, 8157, 8185, 8213, 8241, 8269, 8300, 8325, 8353, 8381, 8409, 8437, 8465, 8493, 8521, 8549, 8577, 8605, 8641, 8669, 8697, 8725, 8753, 8781, 8809, 8837, 8865, 8893, 8921, 8949, 8977, 9005, 9033, 9061, 9089, 9117, 9145, 9173, 9201, 9229, 9257, 9285, 9313, 9344, 9369, 9397, 9425, 9453, 9481, 9509, 9537, 9565, 9593, 9621]
        scanData[Dim.X] = [36.45, 34.05, 31.95, 29.55, 26.85, 24.6, 20.1, 15.6, 9.3, 3.75, -0.45, -6.45, -13.5, -24.15, -31.5, -38.1, -44.55, -51.15, -59.1, -65.1, -71.4, -76.65, -82.05, -85.5, -87.3, -87.9, -87.3, -87.6, -85.5, -83.85, -78.3, -74.25, -70.05, -66.6, -62.85, -58.2, -53.25, -46.35, -39.3, -33, -24.75, -15.9, -5.55, 2.55, 9, 13.8, 17.85, 21.75, 25.2, 28.8, 31.65, 34.2, 36.75, 37.5, 35.55, 30.3, 26.4, 22.2, 21, 17.85, 13.35, 5.85, -1.2, -8.85, -16.8, -25.8, -33.9, -41.4, -48.9, -55.65, -62.55, -68.1, -71.85, -75.9, -79.05, -82.95, -85.8, -87.45, -86.25, -85.05, -82.05, -80.1, -77.1, -73.95, -71.7, -66.6, -61.95, -53.85, -47.25, -41.1, -33.3, -26.4, -16.95, -10.35, -4.35, 0, 4.5, 12.45, 18.9, 26.1, 30, 32.55, 34.05, 34.2, 34.95, 32.4, 31.8, 29.4, 28.5, 25.05, 18.75, 11.7, 4.8, -0.6, -7.8, -16.8, -25.05, -32.7, -38.1, -43.5, -49.35, -55.5, -62.85, -69.3, -75.6, -80.7, -83.7, -85.65, -86.25, -86.7, -87.15, -85.95, -83.1, -79.95, -75.75, -70.5, -67.8, -62.7, -60.45, -51.75, -45, -37.05, -28.2]
        scanData[Dim.Y] = [-387.75, -389.55, -385.2, -385.95, -386.1, -384.15, -381.15, -379.8, -381, -380.1, -380.7, -380.85, -382.95, -384.45, -385.5, -383.85, -381.15, -377.55, -376.2, -375.45, -375, -374.25, -372.6, -373.35, -373.05, -373.35, -373.35, -374.7, -376.5, -378.3, -378.15, -377.7, -375.9, -375.75, -375.45, -376.35, -377.55, -376.95, -376.95, -376.5, -378.6, -379.2, -381.3, -381.45, -382.95, -382.05, -382.2, -381.15, -381.45, -383.1, -384.9, -386.7, -385.5, -384.6, -382.8, -382.2, -382.8, -384.15, -385.5, -385.35, -387.3, -385.35, -382.8, -379.65, -380.7, -382.35, -381.45, -380.7, -380.85, -380.25, -379.35, -377.7, -377.85, -376.95, -376.95, -376.05, -376.8, -375.9, -376.8, -376.8, -377.85, -378.9, -379.8, -379.5, -377.85, -375.9, -376.95, -377.7, -380.55, -378.9, -378, -378.15, -379.35, -382.05, -380.55, -380.4, -381, -383.7, -387.15, -387.15, -387, -385.05, -384.15, -384.3, -384.3, -384.6, -383.25, -384, -384.45, -385.8, -384.6, -384, -382.35, -381.6, -380.1, -380.1, -381.3, -384, -382.8, -382.35, -379.5, -379.2, -378.9, -378.9, -379.05, -378.9, -377.1, -376.65, -377.7, -378.9, -377.1, -374.1, -373.05, -374.55, -374.1, -375.75, -377.1, -378.9, -378, -377.4, -376.8, -376.65, -375.75]
        scanData[Dim.Z] = [-145.05, -149.7, -155.25, -157.8, -159.45, -162.9, -165.9, -171.15, -173.4, -178.8, -179.4, -181.2, -181.8, -183.3, -184.35, -184.35, -182.1, -180, -175.5, -170.85, -164.7, -160.95, -156.45, -150.15, -144.6, -138.3, -132.3, -124.5, -115.95, -109.65, -99.9, -93.6, -84.3, -78.9, -75.15, -72.6, -68.25, -64.05, -63.6, -63.45, -64.35, -63.15, -64.8, -67.35, -72.45, -79.5, -84, -89.4, -95.55, -103.65, -110.55, -115.95, -122.1, -130.2, -139.05, -146.25, -151.95, -156.75, -163.05, -169.05, -171.9, -173.4, -177.3, -182.7, -184.35, -183.45, -181.05, -181.65, -178.05, -176.7, -172.95, -170.1, -165.45, -159.75, -154.35, -147, -140.1, -134.55, -127.95, -121.95, -111.75, -103.05, -97.95, -92.1, -88.35, -81.3, -75.9, -72, -67.65, -68.55, -65.55, -64.65, -63, -63.9, -67.8, -72.6, -76.65, -80.55, -85.05, -94.95, -102.9, -109.8, -116.25, -122.7, -130.05, -136.5, -143.85, -150.6, -156, -162.15, -168.3, -172.35, -175.65, -179.25, -180.6, -184.2, -183.15, -184.95, -183.15, -181.2, -178.95, -174.6, -169.65, -166.2, -162.15, -156.75, -150.75, -143.7, -138.9, -132.6, -125.4, -118.35, -108, -102.45, -94.2, -87.75, -79.35, -73.05, -71.25, -68.55, -66.3, -63.75, -63.9]
    }

    // EXPORTED USER INTERFACES   

    /** 
     * Assuming the buggy is currently spinning clockwise on the spot, capture a 
     * time-stamped sequence of magnetometer readings, from which to set up the compass.
     * NOTE that once scanning is complete, the heading.setNorth() function must then 
     * be called (to process the scanned data) before heading.degrees() will work.
     * @param ms scanning-time in millisecs (long enough for more than one full rotation) 
    */
    //% block="scan for (ms) $ms" 
    //% inlineInputMode=inline 
    //% ms.shadow="timePicker" 
    //% ms.defl=0 
    //% weight=90 

    export function scan(ms: number) {
        // Magnetometer readings are scanned into four internal arrays: times[], xVals[], yVals[] & zVals[],
        // sampled every ~30 ms over the specified duration (usually at least a second)
        if (testing) {
            simulateScan()
            basic.pause(ms)
            return
        }

        // (To cure jitter, each reading is always a rolling sum of SEVEN consecutive readings!)
        let now = input.runningTime()
        let finish = now + ms
        let sum = 0
        let xRoll: number[] = []
        let yRoll: number[] = []
        let zRoll: number[] = []
        // take the first six readings...
        for (let i = 0; i<6; i++) {
            xRoll.push(input.magneticForce(0))
            yRoll.push(input.magneticForce(1))
            zRoll.push(input.magneticForce(2))
            basic.pause(25)
        }

        // continue cranking out rolling sums, adding a new reading and dropping the oldest
        while (now < finish) {
            now = input.runningTime()
            scanTimes.push(now-100) // the time of the middle readings (roughly)
            basic.pause(25)
            xRoll.push(input.magneticForce(0))
            yRoll.push(input.magneticForce(1))
            zRoll.push(input.magneticForce(2))
            sum = 0
            xRoll.forEach(a => sum += a)
            xRoll.shift()
            scanData[Dim.X].push(sum)
            sum = 0
            yRoll.forEach(a => sum += a)
            yRoll.shift()
            scanData[Dim.Y].push(sum)
            sum = 0
            zRoll.forEach(a => sum += a)
            zRoll.shift()
            scanData[Dim.Z].push(sum)
        }
    }


    /**
     * Analyse scanned data to prepare for reading compass-headings.
     * Returns either the spin RPM, or a negative error code:
     * 
     *      -1 : NOT ENOUGH SCAN DATA
     * 
     *      -2 : NOT ENOUGH SCAN ROTATION
     * 
     *      -3 : FIELD STRENGTH TOO WEAK
     */
    //% block="set North" 
    //% inlineInputMode=inline 
    //% weight=80 
    export function setNorth(): number {
    /* 
    As the buggy spins, the magnetic field-vector sweeps out a cone. In the fully general case, 
    this projects onto the plane of each pair of orthogonal axes (XY,YZ,ZX) as an ellipse
    with a certain eccentricity. We will get the best heading discrimination from the plane 
    with the least eccentric ellipse, and having selected those two axes, we'll need to 
    transform readings around the ellipse so that they lie on a circle, giving a relative 
    angle that can (eventually) be offset by a fixed bias to return the true heading.
    */

        // we need at least a second's worth of scanned readings...
        if (scanTimes.length < 40) {
            return -1 // "NOT ENOUGH SCAN DATA"
        }
        return 0
        // find the raw extrema
        let xlo = 999
        let ylo = 999
        let zlo = 999
        let xhi = -999
        let yhi = -999
        let zhi = -999
        let v = 0
        for (let i = 0; i <= scanTimes.length; i++) {
            v = scanData[Dim.X][i]
            if (v < xlo) xlo = v
            if (v > xhi) xhi = v
            v = scanData[Dim.Y][i]
            if (v < ylo) ylo = v
            if (v > yhi) yhi = v
            v = scanData[Dim.Z][i]
            if (v < zlo) zlo = v
            if (v > zhi) zhi = v
        }
        let xOff = (xhi + xlo) / 2
        let yOff = (yhi + ylo) / 2
        let zOff = (zhi + zlo) / 2
        let xAmp = (xhi - xlo) / 2
        let yAmp = (yhi - ylo) / 2
        let zAmp = (zhi - zlo) / 2

        if ((xAmp + yAmp + zAmp) < MarginalField) {
            return -3  // "FIELD STRENGTH TOO WEAK"
        }

        // now find the extreme radii for each axis-pair 
        let xylo = 99999
        let yzlo = 99999
        let zxlo = 99999
        let xyhi = -99999
        let yzhi = -99999
        let zxhi = -99999
        let xyMax: number[] = [0]
        let yzMax: number[] = [0]
        let zxMax: number[] = [0]
        let period = 0
        let x = 0
        let y = 0
        let z = 0
        let xsq = 0
        let ysq = 0
        let zsq = 0
        let xya = 0
        let yza = 0
        let zxa = 0

        for (let i = 0; i <= scanTimes.length; i++) {
            x = scanData[Dim.X][i] - xOff
            y = scanData[Dim.Y][i] - yOff
            z = scanData[Dim.Z][i] - zOff
            xsq = x * x
            ysq = y * y
            zsq = z * z
            v = xsq + ysq // radius-squared is sum of squares
            if (v < xylo) xylo = v
            if (v > xyhi) {
                xyhi = v
                xya = Math.atan2(x, y)
                // record time if new extreme radius found at least ~15 readings on from last one
                let last = scanTimes[xyMax[xyMax.length - 1]]
                if (scanTimes[i] - last > 300) xyMax.push(i)
            }
            v = ysq + zsq
            if (v < yzlo) yzlo = v
            if (v > yzhi) {
                yzhi = v
                yza = Math.atan2(y, z) 
                let last = scanTimes[yzMax[yzMax.length - 1]]
                if (scanTimes[i] - last > 300) yzMax.push(i)
            }

            v = zsq + xsq
            if (v < zxlo) zlo = v
            if (v > zxhi) {
                zxhi = v
                zxa = Math.atan2(z, x) 
                let last = scanTimes[zxMax[zxMax.length - 1]]
                if (scanTimes[i] - last > 300) zxMax.push(i)
            }
        }
        
        if (axes[vDim].amp < MarginalField) {
            return -3  // "FIELD STRENGTH TOO WEAK"
        }
        // use eccentricities to select best axes to use
        let xye = Math.sqrt(xyhi / xylo)
        let yze = Math.sqrt(yzhi / yzlo)
        let zxe = Math.sqrt(zxhi / zxlo)
        if (xye < yze) { // not YZ
            if (xye < zxe) { // not ZX either, so use XY
                plane = "XY"
                uDim = Dim.X
                vDim = Dim.Y
                uOff = xOff
                vOff = yOff
                theta = xya
                scale = xye
                period = 2 * (xyMax.pop() - xyMax[0]) / xyMax.length
            }
        } else { // not XY: either YZ or ZX
            if (yze < zxe) { // not ZX so use YZ
                plane = "YZ"
                uDim = Dim.Y
                vDim = Dim.Z
                uOff = yOff
                vOff = zOff
                theta = yza
                scale = yze
                period = 2 * (yzMax.pop() - yzMax[0]) / yzMax.length
            } else { // not YZ so use ZX
                plane = "ZX"
                uDim = Dim.Z
                vDim = Dim.X
                uOff = zOff
                vOff = xOff
                theta = zxa
                scale = zxe
                period = 2 * (zxMax
    .pop() - zxMax
    [0]) / zxMax
    .length
            }
        }

// We have set up the projection parameters. Now we need to relate them to North.
// Take the average of three new readings to get a stable fix on the current heading
        let uRaw = 0
        let vRaw = 0
        if (testing) { // choose some arbitrary reading for North (we know U=X & V=Z)
            uRaw = scanData[Dim.X][40]
            vRaw = scanData[Dim.Z][40]
        } else {
            // get a rolling sum of three readings
            uRaw = input.magneticForce(uDim)
            vRaw = input.magneticForce(vDim)
            basic.pause(5)
            uRaw += input.magneticForce(uDim)
            vRaw += input.magneticForce(vDim)
            basic.pause(5)
            uRaw += input.magneticForce(uDim)
            vRaw += input.magneticForce(vDim)
        }
        toNorth = project(uRaw/3, vRaw/3)


        // For a clockwise scan, the maths requires the U-dim to lead the V-dim by 90 degrees
        // From the point of view of a buggy spinning clockwise from ~NW, the North vector appears 
        // to rotate anticlockwise, passing the +V axis first, and then the -U axis.
        // Check the timings of their first limits and, if necessary, swap the major/minor dimensions:
        /* if (axes[uDim].time0 < axes[vDim].time0) {
            let temp = uDim
            uDim = vDim
            vDim = temp
        }

        // Also check the polarities of these first limits in case the microbit
        // is mounted backwards: we expect the first uVal<0 and the first vVal>0
        uFlip = -(axes[uDim].limit0 / Math.abs(axes[uDim].limit0)) // = -1 if uVal>0
        vFlip = axes[vDim].limit0 / Math.abs(axes[vDim].limit0)    // = -1 if vVal<0
*/
        // set up datalogger
        datalogger.setColumnTitles("uRaw", "vRaw", "u", "v", "uNew", "vNew", "vScaled", "a")
        datalogger.includeTimestamp(FlashLogTimeStampFormat.None)
        // return average RPM of original scan    
        return 60000 / period
 
    }
    
    // correct a two-axis raw magnetometer reading from the off-centre projected ellipse
    // back onto a centred circle of headings
    function project(uRaw: number, vRaw: number) {
        // recentre point
        let u = uRaw - uOff
        let v = vRaw - vOff
        // rotate by the major-axis angle theta (check direction!)
        let uNew = u * Math.cos(theta) - v * Math.sin(theta)
        let vNew = u * Math.sin(theta) + v * Math.cos(theta)
        // scale V to match U
        let vScaled = vNew * scale
        // return projected angle
        let a = Math.atan2(uNew,vScaled)
        
        datalogger.log(datalogger.createCV("dim", uDim),
            datalogger.createCV("bias", axes[uDim].bias),
            datalogger.createCV("amp", axes[uDim].amp),
            datalogger.createCV("flip", uFlip))    
    }



        
        



        if (axes[vDim].amp < MarginalField) {
            return -3  // "FIELD STRENGTH TOO WEAK"
        } else {
        // return best average RPM of original scan    
            return 120000 / (axes[uDim].period + axes[vDim].period)
        }
    }
  


    /**
     * Read the magnetometer and return the current heading of the buggy in degrees
     */
    //% block="heading" 
    //% inlineInputMode=inline 
    //% weight=70
    export function degrees(): number {
        // read the magnetometer (three times) and return the current heading in degrees from North
        let uRaw = 0
        let vRaw = 0
        if (testing) { // NOTE: a-priori knowledge U=X & V=Z !
            uRaw = scanData[Dim.X][test]
            vRaw = scanData[Dim.Z][test]
            test += 4
            if (test > scanData[Dim.X].length - 2) test = 0 // roll round
        } else {
            // get a rolling sum of three readings
            uRaw = input.magneticForce(uDim)
            vRaw = input.magneticForce(vDim)
            basic.pause(5)
            uRaw += input.magneticForce(uDim)
            vRaw += input.magneticForce(vDim)
            basic.pause(5)
            uRaw += input.magneticForce(uDim)
            vRaw += input.magneticForce(vDim)
        }
        // project reading from ellipse to circle and relate to North
        // (converting from radians to degrees)
        let val = 57.29578 * (project(uRaw,vRaw) - toNorth)
        // shift negative [-180...0] range to positive [180...360]
        val = (val + 360) % 360
        datalogger.log(datalogger.createCV("uRaw", uRaw),
            datalogger.createCV("vRaw", vRaw),
            datalogger.createCV("val", val))
        return val
    }


    //% block="switch test mode" 
    //% inlineInputMode=inline 
    //% weight=50
    export function testMode(turnOn: boolean) {
        testing = turnOn
    }
}
