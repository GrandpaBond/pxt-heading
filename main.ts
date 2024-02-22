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
        nLimits: number // tue count of local extremes detected
        time0: number   // the timestamp of the first limit == Limit[0].time
        peak0: number   // the value of the first limit == Limit[0].value
        period: number  // the RPM of the scanned wave-form in this axis

       constructor(dim: Dim) {
       // dim should be redundant, as it always matches index as axes[] array is built
           this.dim = dim
           this.limits = []
        }

        // Method to characterise the scanData for each axis
        characterise(timeStamp: number[], scanData: number[]) {
        // First we extract local maxima and minima from a sequence of scanned wave-form readings
        // For element [d] in [...a,b,c,d,e...], we multiply the averaged slope behind it
        // (given by [d]-[a]), with the slope ahead (given by [e]-[b]). 
        // Whenever the sign of the slopechanges (i.e this product is negative) we record a Limit: 
        // a duple object comprising the peak value [d], together with its time-stamp.
            this.nLimits = 0
            let then = timeStamp[0] - 201
            // loop through from [d] onwards
            for (let i = 3; i <= timeStamp.length - 3; i++) {
                let ahead = (scanData[i+2] - scanData[i-1])
                let behind = (scanData[i+1] - scanData[i-2])
            // An inflection-point is where the product will be negative.
            // (Ignore any too-close crossings, arising from excessive jitter)
                if (((ahead * behind) <= 0) 
                && ((timeStamp[i] - then) > 200)) {
                    let newLimit = new Limit
                    newLimit.value = scanData[i]
                    newLimit.time = timeStamp[i]
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

        // work out the periodicity
            let spans = this.nLimits - 1
            this.period = 2* (this.limits[spans].time - this.limits[0].time) / spans
        }
    }

    // GLOBALS
    const MarginalField = 50 // minimum acceptable amplitude for magnetometer readings

    let scanTimes: number[] = [] // sequence of time-stamps for scanned readings 
    let xScanData: number[] = [] // scanned sequence of magnetometer X-Axis readings 
    let yScanData: number[] = [] // scanned sequence of magnetometer Y-Axis readings
    let zScanData: number[] = [] // scanned sequence of magnetometer Z-Axis readings  
    
    let axes: Axis[] = []
    axes.push(new Axis(Dim.X))
    axes.push(new Axis(Dim.Y))
    axes.push(new Axis(Dim.Z))
    
    let uDim = -1 // will hold the horizontal axis for North vector
    let vDim = -1 // will hold the vertical axis for North vector
    
    let uFlip = 1 // set to -1 if uDim polarity is inverted
    let vFlip = 1 // set to -1 if vDim polarity is inverted
    let plane = "??"



    // USER INTERFACE

    // use some sample data, while debugging...
    export function simulateScan() {    
        scanTimes = [5637, 5665, 5693, 5721, 5749, 5777, 5805, 5833, 5861, 5889, 5917, 5945, 5973, 6001, 6029, 6057, 6087, 6113, 6141, 6169, 6197, 6225, 6253, 6281, 6309, 6337, 6365, 6393, 6421, 6449, 6477, 6505, 6533, 6561, 6589, 6617, 6645, 6673, 6701, 6729, 6757, 6785, 6813, 6841, 6869, 6897, 6925, 6956, 6981, 7009, 7037, 7065, 7093, 7121, 7149, 7177, 7205, 7233, 7261, 7289, 7317, 7345, 7373, 7401, 7429, 7460, 7485, 7513, 7541, 7569, 7597, 7625, 7653, 7681, 7709, 7737, 7765, 7793, 7821, 7849, 7877, 7905, 7933, 7961, 7989, 8017, 8045, 8073, 8101, 8129, 8157, 8185, 8213, 8241, 8269, 8300, 8325, 8353, 8381, 8409, 8437, 8465, 8493, 8521, 8549, 8577, 8605, 8641, 8669, 8697, 8725, 8753, 8781, 8809, 8837, 8865, 8893, 8921, 8949, 8977, 9005, 9033, 9061, 9089, 9117, 9145, 9173, 9201, 9229, 9257, 9285, 9313, 9344, 9369, 9397, 9425, 9453, 9481, 9509, 9537, 9565, 9593, 9621]
        xScanData = [36.45, 34.05, 31.95, 29.55, 26.85, 24.6, 20.1, 15.6, 9.3, 3.75, -0.45, -6.45, -13.5, -24.15, -31.5, -38.1, -44.55, -51.15, -59.1, -65.1, -71.4, -76.65, -82.05, -85.5, -87.3, -87.9, -87.3, -87.6, -85.5, -83.85, -78.3, -74.25, -70.05, -66.6, -62.85, -58.2, -53.25, -46.35, -39.3, -33, -24.75, -15.9, -5.55, 2.55, 9, 13.8, 17.85, 21.75, 25.2, 28.8, 31.65, 34.2, 36.75, 37.5, 35.55, 30.3, 26.4, 22.2, 21, 17.85, 13.35, 5.85, -1.2, -8.85, -16.8, -25.8, -33.9, -41.4, -48.9, -55.65, -62.55, -68.1, -71.85, -75.9, -79.05, -82.95, -85.8, -87.45, -86.25, -85.05, -82.05, -80.1, -77.1, -73.95, -71.7, -66.6, -61.95, -53.85, -47.25, -41.1, -33.3, -26.4, -16.95, -10.35, -4.35, 0, 4.5, 12.45, 18.9, 26.1, 30, 32.55, 34.05, 34.2, 34.95, 32.4, 31.8, 29.4, 28.5, 25.05, 18.75, 11.7, 4.8, -0.6, -7.8, -16.8, -25.05, -32.7, -38.1, -43.5, -49.35, -55.5, -62.85, -69.3, -75.6, -80.7, -83.7, -85.65, -86.25, -86.7, -87.15, -85.95, -83.1, -79.95, -75.75, -70.5, -67.8, -62.7, -60.45, -51.75, -45, -37.05, -28.2]
        yScanData = [-387.75, -389.55, -385.2, -385.95, -386.1, -384.15, -381.15, -379.8, -381, -380.1, -380.7, -380.85, -382.95, -384.45, -385.5, -383.85, -381.15, -377.55, -376.2, -375.45, -375, -374.25, -372.6, -373.35, -373.05, -373.35, -373.35, -374.7, -376.5, -378.3, -378.15, -377.7, -375.9, -375.75, -375.45, -376.35, -377.55, -376.95, -376.95, -376.5, -378.6, -379.2, -381.3, -381.45, -382.95, -382.05, -382.2, -381.15, -381.45, -383.1, -384.9, -386.7, -385.5, -384.6, -382.8, -382.2, -382.8, -384.15, -385.5, -385.35, -387.3, -385.35, -382.8, -379.65, -380.7, -382.35, -381.45, -380.7, -380.85, -380.25, -379.35, -377.7, -377.85, -376.95, -376.95, -376.05, -376.8, -375.9, -376.8, -376.8, -377.85, -378.9, -379.8, -379.5, -377.85, -375.9, -376.95, -377.7, -380.55, -378.9, -378, -378.15, -379.35, -382.05, -380.55, -380.4, -381, -383.7, -387.15, -387.15, -387, -385.05, -384.15, -384.3, -384.3, -384.6, -383.25, -384, -384.45, -385.8, -384.6, -384, -382.35, -381.6, -380.1, -380.1, -381.3, -384, -382.8, -382.35, -379.5, -379.2, -378.9, -378.9, -379.05, -378.9, -377.1, -376.65, -377.7, -378.9, -377.1, -374.1, -373.05, -374.55, -374.1, -375.75, -377.1, -378.9, -378, -377.4, -376.8, -376.65, -375.75]
        zScanData = [-145.05, -149.7, -155.25, -157.8, -159.45, -162.9, -165.9, -171.15, -173.4, -178.8, -179.4, -181.2, -181.8, -183.3, -184.35, -184.35, -182.1, -180, -175.5, -170.85, -164.7, -160.95, -156.45, -150.15, -144.6, -138.3, -132.3, -124.5, -115.95, -109.65, -99.9, -93.6, -84.3, -78.9, -75.15, -72.6, -68.25, -64.05, -63.6, -63.45, -64.35, -63.15, -64.8, -67.35, -72.45, -79.5, -84, -89.4, -95.55, -103.65, -110.55, -115.95, -122.1, -130.2, -139.05, -146.25, -151.95, -156.75, -163.05, -169.05, -171.9, -173.4, -177.3, -182.7, -184.35, -183.45, -181.05, -181.65, -178.05, -176.7, -172.95, -170.1, -165.45, -159.75, -154.35, -147, -140.1, -134.55, -127.95, -121.95, -111.75, -103.05, -97.95, -92.1, -88.35, -81.3, -75.9, -72, -67.65, -68.55, -65.55, -64.65, -63, -63.9, -67.8, -72.6, -76.65, -80.55, -85.05, -94.95, -102.9, -109.8, -116.25, -122.7, -130.05, -136.5, -143.85, -150.6, -156, -162.15, -168.3, -172.35, -175.65, -179.25, -180.6, -184.2, -183.15, -184.95, -183.15, -181.2, -178.95, -174.6, -169.65, -166.2, -162.15, -156.75, -150.75, -143.7, -138.9, -132.6, -125.4, -118.35, -108, -102.45, -94.2, -87.75, -79.35, -73.05, -71.25, -68.55, -66.3, -63.75, -63.9]
    }

    // EXPORTED USER INTERFACES   
  
    /** 
     * Assuming the buggy is currently spinning clockwise on the spot, capture a 
     * time-stamped sequence of magnetometer readings, from which to set up the compass.
     * NOTE that once scanning is complete, the heading.analyse() function must then 
     * be called (to process the scanned data) before heading.degrees() will work.
     * The buggy should be pointing approximately NW to start with, so that it passes
     * the cardinal points in the order {N, E, S, W}.
     * @param ms scanning-time in millisecs (long enough for at least one full rotation) 
    */ 
     //% block="scan for (ms) $ms" 
     //% inlineInputMode=inline 
     //% ms.shadow="timePicker" 
     //% ms.defl=0 
     //% weight=90 
     
    export function scan(ms: number) {
    // Magnetometer readings are scanned into four internal arrays: times[], xVals[], yVals[] & zVals[],
    // sampled every ~30 ms over the specified duration (usually at least a second)
    // (To cure jitter, each reading is always a rolling sum of three consecutive readings)
        let now = input.runningTime()
        let finish = now + ms
        let sum = 0
        let xRoll: number[] = []
        let yRoll: number[] = []
        let zRoll: number[] = []
        // take the first couple of readings...
        xRoll.push(input.magneticForce(0))
        yRoll.push(input.magneticForce(1))
        zRoll.push(input.magneticForce(2))
        basic.pause(25)
        xRoll.push(input.magneticForce(0))
        yRoll.push(input.magneticForce(1))
        zRoll.push(input.magneticForce(2))
        // continue cranking out rolling sums, adding a new reading and dropping the oldest
        while (now < finish) {
            now = input.runningTime()
            scanTimes.push(now) // the time of the middle readings (roughly)
            basic.pause(25)
            xRoll.push(input.magneticForce(0))
            yRoll.push(input.magneticForce(1))
            zRoll.push(input.magneticForce(2))
            sum = 0
            xRoll.forEach(a => sum += a)
            xRoll.shift() 
            xScanData.push(sum)
            sum = 0
            yRoll.forEach(a => sum += a)
            yRoll.shift()
            yScanData.push(sum)
            sum = 0
            zRoll.forEach(a => sum += a)
            zRoll.shift()   
            zScanData.push(sum)  
        }
    }


    /**
     * Analyse scanned data to prepare for reading compass-headings.
     * returns either the spin RPM, or a negative error code:
     *      -1 : NOT ENOUGH SCAN DATA
     *      -2 : NOT ENOUGH SCAN ROTATION
     *      -3 : FIELD STRENGTH TOO WEAK
     */
    //% block="analyse scan" 
    //% inlineInputMode=inline 
    //% weight=80 
    export function analyseScan(): number {
    // Analyse the arrays of scanned data to:
    // a) assign major and minor magnetometer dimensions [u,v] and find North
    // b) calculate the normalisation offsets & scaling to be applied to [u,v]
    // c) detect the periodicity
    // returns the (quite interesting) RPM of the original scan (or a negative error-code)

        // we need at least a second's worth of readings...
        if (scanTimes.length < 40) {
            return -1 // "NOT ENOUGH SCAN DATA"
        }

        // datalogger.includeTimestamp(FlashLogTimeStampFormat.None)
        // datalogger.mirrorToSerial(true)
        // datalogger.setColumnTitles("xt", "xlim", "yt", "ylim", "zt", "zlim", )

        axes[0].characterise(scanTimes, xScanData)
        axes[1].characterise(scanTimes, yScanData)
        axes[2].characterise(scanTimes, zScanData)
        

        // To guarantee an example of a maximum and a minimum in each dimension, we
        // need to have scanned at least one and a bit full revolutions of the buggy...
        if ((axes[0].nLimits < 3) || (axes[1].nLimits < 3) || (axes[2].nLimits < 3)) {
            return -2  // "NOT ENOUGH SCAN ROTATION"
        }

        // Choose the two dimensions showing the biggest amplitudes. From now on, we work
        // with the projection of the magnetic flux vector onto the selected [U,V] plane.
        let hi = Math.max(axes[Dim.X].amp, Math.max(axes[Dim.Y].amp, axes[Dim.Z].amp))
        let lo = Math.min(axes[Dim.X].amp, Math.min(axes[Dim.Y].amp, axes[Dim.Z].amp))
        uDim = axes.find(i => i.amp === hi).dim
        vDim = axes.find(i => (i.amp != hi) && (i.amp != lo)).dim

        if (axes[vDim].amp < MarginalField) {
            return -3  // "FIELD STRENGTH TOO WEAK"
        }

        // For a clockwise scan, the maths requires the U-dim to lead the V-dim by 90 degrees
        // (so that e.g. when U = 2*V and both are positive we get +60 degrees).
        // From the point of view of a buggy spinning clockwise from ~NW, the North vector appears 
        // to rotate anticlockwise, passing the +V axis first, and then the -U axis.
        // Check the timings of their first limits and, if necessary, swap the major/minor dimensions:
        if (axes[uDim].limits[0].time < axes[vDim].limits[0].time) {
            let temp = uDim
            uDim = vDim
            vDim = temp
        }
        let uVal = axes[uDim].limits[0].value
        let vVal = axes[vDim].limits[0].value
        // Also check the polarities of these first limits in case the microbit 
        // is mounted backwards: we expect uVal<0 and vVal>0
        uFlip = -(uVal / Math.abs(uVal)) // = 1 unless uVal>0
        vFlip = vVal / Math.abs(vVal)    // = 1 unless vVal<0

        // return best average RPM of original scan    
        return 120000 / (axes[uDim].period + axes[vDim].period)
    }
    /**
     * Read the magnetometer.
     * returns the current heading of the buggy in degrees
     */
    //% block="heading" 
    //% inlineInputMode=inline 
    //% weight=70
    export function degrees(): number {
    // read the magnetometer (three times) and return the current heading in degrees from North
        let uRaw = 0
        let vRaw = 0
        if (testing) { // NOTE: a-priori knowledge U=X & V=Z !
            uRaw = xScanData[test]
            vRaw = zScanData[test]
            test += 4
            if (test > xScanData.length-2) test = 0 // roll round
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
        // normalise the horizontal & vertical components
        let uPart = uFlip * (uRaw - axes[uDim].bias) / axes[uDim].amp
        let vPart = vFlip * (vRaw - axes[vDim].bias) / axes[vDim].amp
        let val = 57.29578 * Math.atan2(uPart, vPart)
        // shift negative [-180...0] range to positive [180...360]
        val = (val + 360) % 360
        datalogger.log(datalogger.createCV("uRaw", uRaw),
            datalogger.createCV("vRaw", vRaw),
            datalogger.createCV("u", uPart),
            datalogger.createCV("v", vPart),
            datalogger.createCV("val", val))
        return val
    }
    
    export function dumpData() {
        datalogger.deleteLog()
        datalogger.setColumnTitles("t","x","y","z")
        for (let i = 0; i < scanTimes.length; i++) {
            datalogger.log(datalogger.createCV("t", scanTimes[i]),
                            datalogger.createCV("x", xScanData[i]),
                            datalogger.createCV("y", yScanData[i]),
                            datalogger.createCV("z", zScanData[i]))
        }
    }
}
