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
        analyse(timeStamp: number[], scanData: number[]) {
        // First we extract local maxima and minima from a sequence of scanned wave-form readings
        // For element [d] in [...a,b,c,d,e,f,g...], we multiply the averaged slope behind it
        // (given by [d]-[a]), with the slope ahead (given by [g]-[d]). 
        // Whenever the sign of the slopechanges (i.e this product is negative) we record a Limit: 
        // a duple object comprising the peak value [d], together with its time-stamp.
            this.nLimits = 0
            let then = timeStamp[0] - 100
            // loop through from [d] onwards
            for (let i = 3; i <= timeStamp.length - 4; i++) {
            // An inflection-point is where the product will be negative.
            // (Ignore any too-close crossings, arising from excessive jitter)
                if ((scanData[i] - scanData[i-3] * (scanData[i+3] - scanData[i]) < 0) && ((timeStamp[i] - then) > 100)) {
                    let newLimit = new Limit
                    newLimit.value = scanData[i]
                    newLimit.time = timeStamp[i]
                    then = timeStamp[i]
                    this.limits.push(newLimit)
                    this.nLimits++
                }
            }

        // Use the averages of paired limits to set the magnitude and offset.
        // Needs at least two of each (resulting from a complete spin)
            this.amp = 0
            this.bias = 0
            if (this.nLimits > 2) {
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
    
    let quadOffset = 0
    let reversed = false
    let plane = "??"



    // USER INTERFACE

    // use some sample data, while debugging...
    export function simulateScan() {
        scanTimes = [6193, 6229, 6257, 6285, 6313, 6341, 6369, 6397, 6425, 6453, 6481, 6509, 6537, 6565, 6593, 6621, 6649, 6677, 6705, 6733, 6761, 6789, 6817, 6845, 6873, 6901, 6929, 6957, 6985, 7013, 7041, 7069, 7097, 7125, 7153, 7181, 7209, 7237, 7265, 7293, 7321, 7349, 7377, 7405, 7433, 7461, 7489, 7517, 7545, 7573, 7601, 7629, 7657, 7685, 7713, 7741, 7769, 7797, 7825, 7853, 7881, 7909, 7937, 7965, 7993, 8021, 8049, 8077, 8105, 8133, 8161, 8189, 8217, 8245, 8273, 8301, 8329, 8357, 8385, 8413, 8441, 8469, 8497, 8525, 8553, 8581, 8609, 8637, 8665, 8693, 8721, 8749, 8777, 8805, 8833, 8861, 8889, 8917, 8945, 8973, 9001, 9029, 9057, 9085, 9113, 9141, 9169, 9197, 9225, 9253, 9281, 9309, 9337, 9365, 9393, 9421, 9449, 9477, 9505, 9533, 9561, 9589, 9617, 9645, 9673, 9701, 9729, 9757, 9785, 9813, 9841, 9869, 9897, 9925, 9953, 9981, 10009, 10037, 10065]
        xScanData = [-6.75, -5.4, -4.5, 1.05, 2.25, 5.1, 7.65, 10.8, 13.5, 15.9, 15.6, 16.05, 17.7, 18.15, 18.75, 19.35, 18.6, 18.15, 15.9, 16.65, 13.95, 13.05, 10.5, 8.55, 5.85, 2.4, 0.6, -3, -4.35, -8.1, -8.85, -9.9, -12.6, -14.7, -18.6, -19.5, -19.65, -22.35, -22.2, -22.8, -22.2, -22.35, -21, -20.7, -20.1, -19.05, -17.25, -16.05, -15.45, -12.15, -9.75, -6.75, -4.65, 0.15, 1.95, 3.3, 5.7, 6.9, 10.5, 12.75, 12.6, 15.3, 15.45, 17.7, 18.15, 19.5, 18.75, 18.75, 16.8, 17.1, 15.6, 14.55, 13.65, 12.6, 12.15, 8.25, 6, 3.3, 0.15, -0.9, -5.25, -7.05, -11.1, -11.7, -14.55, -16.35, -17.55, -18.6, -18.9, -21.75, -22.5, -22.8, -23.85, -24.3, -21.9, -22.35, -20.55, -18.15, -15.75, -13.5, -11.7, -9.3, -7.35, -5.85, -2.1, 1.05, 2.4, 4.65, 7.35, 10.95, 10.95, 13.8, 14.4, 17.85, 17.7, 17.1, 18.9, 18.9, 18.75, 19.35, 17.25, 16.65, 15.6, 14.85, 12.9, 11.55, 9, 6.6, 2.1, 1.05, -3.15, -5.85, -7.5, -10.05, -11.1, -12.75, -15.3, -18, -19.2]
        yScanData = [-126, -127.65, -125.25, -125.4, -127.8, -127.35, -126.9, -127.05, -126.6, -127.05, -126.6, -125.7, -126.9, -126.45, -127.8, -127.65, -127.05, -127.2, -127.2, -126.9, -129.45, -125.25, -126.15, -128.55, -128.25, -127.8, -127.95, -129.9, -126.75, -126.3, -129.45, -128.4, -127.35, -126.6, -128.1, -127.2, -125.85, -126.75, -126.45, -127.35, -126.15, -126.6, -126, -125.85, -126.75, -126.15, -127.05, -126.15, -125.85, -126.45, -126.15, -127.05, -126.3, -126.3, -127.2, -126, -127.65, -124.8, -125.85, -126.45, -127.05, -125.85, -127.5, -129.6, -128.4, -125.4, -128.25, -127.35, -127.05, -127.5, -128.1, -127.05, -127.35, -127.65, -127.65, -127.05, -126.45, -125.4, -128.4, -125.55, -127.2, -126.75, -126, -127.05, -127.5, -126, -130.05, -126.9, -126.75, -126.3, -125.85, -126.15, -127.2, -126.6, -126.15, -125.1, -126.75, -126.75, -126.3, -125.85, -127.5, -127.35, -126.3, -127.2, -126.15, -127.05, -126.6, -126.3, -126.75, -127.8, -127.35, -126.15, -125.7, -126.9, -126.15, -127.65, -128.7, -126.6, -127.65, -128.1, -127.8, -127.35, -127.95, -126.15, -126, -126.9, -126.75, -126.75, -128.1, -127.2, -127.05, -126.3, -125.55, -129.15, -127.2, -127.35, -125.85, -125.85, -129]
        zScanData = [-19.5, -17.4, -18.6, -20.4, -18, -19.65, -21.9, -20.55, -24.15, -25.05, -27.6, -30.15, -29.7, -34.8, -37.35, -39.9, -42.3, -44.85, -47.85, -48, -51.15, -54.75, -55.8, -55.65, -58.65, -56.85, -57.45, -56.55, -59.55, -58.65, -55.2, -57.15, -55.65, -54.6, -52.65, -48.6, -49.8, -46.95, -43.8, -40.65, -41.25, -36.6, -34.05, -32.85, -30, -28.8, -24.75, -23.7, -23.55, -21.75, -19.8, -18.9, -18.9, -18.3, -18.3, -17.85, -19.65, -21.45, -23.85, -23.1, -25.35, -27.15, -27.9, -30.45, -34.95, -39.45, -38.7, -42.3, -45.3, -45.15, -47.1, -49.8, -50.85, -52.5, -54.75, -55.95, -58.8, -58.5, -56.55, -59.1, -59.85, -57.9, -57.75, -55.8, -57.15, -55.35, -51, -49.05, -47.4, -44.25, -42.45, -39.45, -36.6, -36.75, -31.05, -30.45, -29.7, -26.7, -23.7, -23.7, -19.65, -19.35, -19.05, -18.45, -18.15, -18.15, -18.45, -20.25, -20.85, -21.15, -21.75, -24.3, -27.3, -27.6, -31.35, -31.05, -34.8, -37.95, -42.75, -41.7, -45.75, -46.95, -49.2, -53.4, -52.05, -55.2, -56.4, -57.6, -58.65, -57.3, -57.6, -59.1, -58.65, -55.8, -56.7, -53.55, -54, -51, -47.7]
    }


    // Assuming the buggy is currently spinning on the spot, capture a time-stamped sequence
    // of magnetometer readings into four internal arrays: times[], xVals[], yVals[] & zVals[],
    // sampled every ~30 ms over the specified duration (which should be at least a second)
    // A clockwise scan should start with the buggy pointing roughly NW, so that it passes
    // N before E, W or S.
    // (To cure jitter, each reading is always a rolling sum of three consecutive readings)
    
    export function scan(duration: number) {
        let now = input.runningTime()
        let finish = now + duration
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

    // While applying a 5-element rolling average to smooth out jitter,
    // analyse the arrays of scanned data to:
    // a) assign major and minor magnetometer dimensions [u,v] and find North
    // b) calculate the normalisation offsets to be applied to [u,v]
    // c) detect the periodicity and rotation sense
    // returns the (quite interesting) RPM of the original scan (or a negative error-code)
    export function prepare(clockwise: boolean): number {

        // we need at least a second's worth of readings...
        if (scanTimes.length < 40) {
            return -1 // "NOT ENOUGH SCAN DATA"
        }

        // datalogger.includeTimestamp(FlashLogTimeStampFormat.None)
        // datalogger.mirrorToSerial(true)
        // datalogger.setColumnTitles("xt", "xlim", "yt", "ylim", "zt", "zlim", )

        axes[0].analyse(scanTimes, xScanData)
        axes[1].analyse(scanTimes, yScanData)
        axes[2].analyse(scanTimes, zScanData)
        

        // To guarantee an example of a maximum and a minimum in each dimension, we
        // need to hace scanned at least one and a bit full revolutions of the buggy...
        if ((axes[0].nLimits < 2) || (axes[1].nLimits < 2) || (axes[2].nLimits < 2)) {
            return -2  // "NOT ENOUGH SCAN ROTATION"
        }

        // Choose the two dimensions showing the biggest amplitudes. From now on, we work
        // with the projection of the magnetic flux vector onto the selected [U,V] plane.
        let hi = Math.max(axes[Dim.X].amp, Math.max(axes[Dim.Y].amp, axes[Dim.Z].amp))
        let lo = Math.min(axes[Dim.X].amp, Math.max(axes[Dim.Y].amp, axes[Dim.Z].amp))
        uDim = axes.find(i => i.amp === hi).dim
        vDim = axes.find(i => (i.amp != hi) && (i.amp != lo)).dim

        // For a clockwise scan, the maths requires the U-dim to lead the V-dim by 90 degrees
        // (so that e.g. when U = 2*V and both are positive we get +60 degrees).
        // Check the timings of their first limits and, if necessary, swap the major/minor dimensions:
        if (axes[uDim].limits[0].time < axes[vDim].limits[0].time) {
            let temp = uDim
            uDim = vDim
            vDim = temp
        }

        // Also check polarities of these first limits in case the microbit is mounted backwards
        quadOffset = 0
        let uVal = axes[uDim].limits[0].value
        let vVal = axes[vDim].limits[0].value
        if (uVal > vVal) quadOffset += 180
        if (Math.abs(uVal) > Math.abs(vVal)) quadOffset += 90

        // return RPM of original scan    
        return 60000 / axes[uDim].period

    }

    // read the magnetometer (three times) and return the current heading in degrees from North
    export function degrees(): number {
        let uRaw = 0
        let vRaw = 0
        if (testing) { // NOTE: prior knowledge U=X & V=Z !
            uRaw = xScanData[test] + xScanData[test + 1] + xScanData[test+2]
            vRaw = zScanData[test] + zScanData[test + 1] + zScanData[test+2]
            test += 5
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
        let uPart = (uRaw - axes[uDim].bias) / axes[uDim].amp
        let vPart = (vRaw - axes[vDim].bias) / axes[vDim].amp
        let val = 57.29578 * Math.atan2(uPart, vPart) + quadOffset
        if (reversed) {
            val = 360 - val // sensor upside-down
        }
        datalogger.log(datalogger.createCV("uRaw", uRaw),
            datalogger.createCV("vRaw", vRaw),
            datalogger.createCV("u", uPart),
            datalogger.createCV("v", vPart),
            datalogger.createCV("val", val))
        return (val + 360) % 360
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
