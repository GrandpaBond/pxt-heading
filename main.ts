namespace heading {

    // GLOBALS
    let xData: number[] = []
    let yData: number[] = []
    let zData: number[] = []
    let stamp: number[] = []

    let uDim = Dimension.X
    let vDim = Dimension.Z
    let uOff = 0
    let vOff = 0
    let magScale = 1
    let reversed = false
    let plane = "??"

    // a Limit is a [time,value] duple
    class Limit {
        time: number
        value: number
        constructor() {
            this.time = 0
            this.value = 0
        }
    }

    let xLimits: Limit[] = []
    let yLimits: Limit[] = []
    let zLimits: Limit[] = []
    let uLimits: Limit[] = []
    let vLimits: Limit[] = []

    export function simulateScan() {
        // use some sample data, while debugging...
        stamp = [6193, 6229, 6257, 6285, 6313, 6341, 6369, 6397, 6425, 6453, 6481, 6509, 6537, 6565, 6593, 6621, 6649, 6677, 6705, 6733, 6761, 6789, 6817, 6845, 6873, 6901, 6929, 6957, 6985, 7013, 7041, 7069, 7097, 7125, 7153, 7181, 7209, 7237, 7265, 7293, 7321, 7349, 7377, 7405, 7433, 7461, 7489, 7517, 7545, 7573, 7601, 7629, 7657, 7685, 7713, 7741, 7769, 7797, 7825, 7853, 7881, 7909, 7937, 7965, 7993, 8021]
        xData = [-6.75, -5.4, -4.5, 1.05, 2.25, 5.1, 7.65, 10.8, 13.5, 15.9, 15.6, 16.05, 17.7, 18.15, 18.75, 19.35, 18.6, 18.15, 15.9, 16.65, 13.95, 13.05, 10.5, 8.55, 5.85, 2.4, 0.6, -3, -4.35, -8.1, -8.85, -9.9, -12.6, -14.7, -18.6, -19.5, -19.65, -22.35, -22.2, -22.8, -22.2, -22.35, -21, -20.7, -20.1, -19.05, -17.25, -16.05, -15.45, -12.15, -9.75, -6.75, -4.65, 0.15, 1.95, 3.3, 5.7, 6.9, 10.5, 12.75, 12.6, 15.3, 15.45, 17.7, 18.15, 19.5]
        yData = [-126, -127.65, -125.25, -125.4, -127.8, -127.35, -126.9, -127.05, -126.6, -127.05, -126.6, -125.7, -126.9, -126.45, -127.8, -127.65, -127.05, -127.2, -127.2, -126.9, -129.45, -125.25, -126.15, -128.55, -128.25, -127.8, -127.95, -129.9, -126.75, -126.3, -129.45, -128.4, -127.35, -126.6, -128.1, -127.2, -125.85, -126.75, -126.45, -127.35, -126.15, -126.6, -126, -125.85, -126.75, -126.15, -127.05, -126.15, -125.85, -126.45, -126.15, -127.05, -126.3, -126.3, -127.2, -126, -127.65, -124.8, -125.85, -126.45, -127.05, -125.85, -127.5, -129.6, -128.4, -125.4]
        zData = [-19.5, -17.4, -18.6, -20.4, -18, -19.65, -21.9, -20.55, -24.15, -25.05, -27.6, -30.15, -29.7, -34.8, -37.35, -39.9, -42.3, -44.85, -47.85, -48, -51.15, -54.75, -55.8, -55.65, -58.65, -56.85, -57.45, -56.55, -59.55, -58.65, -55.2, -57.15, -55.65, -54.6, -52.65, -48.6, -49.8, -46.95, -43.8, -40.65, -41.25, -36.6, -34.05, -32.85, -30, -28.8, -24.75, -23.7, -23.55, -21.75, -19.8, -18.9, -18.9, -18.3, -18.3, -17.85, -19.65, -21.45, -23.85, -23.1, -25.35, -27.15, -27.9, -30.45, -34.95, -39.45]
    }

    // Use the average max and min values of a list of limits to get the magnitude and offset
    // returns [0,0] unless there are at least two of each (resulting from a complete spin)
    function getRange(limits: Limit[]): [number, number] {
        let amp = 0
        let off = 0
        if (limits.length > 0) {
            // 1st pass to get offset
            let sum = 0
            for (let i = 0; i < limits.length; i++) {
                sum += limits[i].value
            }
            off = sum / limits.length

            // 2nd pass to get amplitude
            sum = 0
            for (let j = 0; j < limits.length; j++) {
                sum += Math.abs(limits[j].value - off)
            }
            amp = sum / limits.length
        }

        return [amp, off]
    }


    // Assuming the buggy is currently spinning on the spot, capture a time-stamped sequence
    // of magnetometer readings into four internal arrays: times[], xVals[], yVals[] & zVals[],
    // sampled every ~30 ms over the specified duration (which should be at least a second)
    export function scan(duration: number) {
        basic.showNumber(xData.length)
        let now = input.runningTime()
        let finish = now + duration
        while (now < finish) {
            xData.push(input.magneticForce(Dimension.X))
            yData.push(input.magneticForce(Dimension.Y))
            zData.push(input.magneticForce(Dimension.Z))
            stamp.push(now)
            basic.pause(25)
            now = input.runningTime()
        }
    }

    // While applying a 5-element rolling average to smooth out jitter,
    // analyse the arrays of scanned data to:
    // a) assign major and minor magnetometer dimensions [u,v] and find North
    // b) calculate the normalisation offsets to be applied to [u,v]
    // c) detect the periodicity and rotation sense
    // returns the (quite interesting) RPM of the original scan (or a negative error-code)
    export function prepare(clockwise: boolean): number {

        let uAmp = 0
        let vAmp = 0
        let uOff2 = 0
        let vOff2 = 0

        // we need at least a second's worth of readings...
        if (stamp.length < 40) {
            return -1 // "NOT ENOUGH SCAN DATA"
        }

        // For each dimension, extract the times and smoothed values of the local maxima and minima.
        // For element c in [...a,b,c,d,e...], we take the sum & difference of (a+b+c) and (c+d+e).
        // Wherever the difference (c+d+e)-(a+b+c) changes sign we record a Limit: a duple object
        // comprising the time-stamp, and the smoothed peak value: (c+d+e)+(a+b+c)

        let diffNow = xData[2]
        let diffWas = diffNow
        let then = stamp[0] - 20

        datalogger.includeTimestamp(FlashLogTimeStampFormat.None)
        datalogger.mirrorToSerial(false)
        datalogger.setColumnTitles("xt", "xlim")

        for (let k = 2; k <= stamp.length - 3; k++) {
            let behind = xData[k] + xData[k - 1] + xData[k - 2]
            let ahead = xData[k] + xData[k + 1] + xData[k + 2]
            diffWas = diffNow
            diffNow = ahead - behind
            // difference crosses zero at an inflection-point, so its sign changes
            // (but ignoring any jittery repeated crossings)
            if ((diffWas * diffNow < 0) && ((stamp[k] - then) > 100)){
                let maybe = new Limit
                maybe.value = ahead + behind // = local average x 6
                maybe.time = stamp[k]
                then = stamp[k]
                xLimits.push(maybe)
                datalogger.log(datalogger.createCV("xt", maybe.time),
                    datalogger.createCV("xlim", maybe.value))
            }
        }

        diffNow = yData[2]
        diffWas = diffNow
        then = stamp[0] - 20
        for (let l = 2; l <= stamp.length - 3; l++) {
            let behind2 = yData[l] + yData[l - 1] + yData[l - 2]
            let ahead2 = yData[l] + yData[l + 1] + yData[l + 2]
            diffWas = diffNow
            diffNow = ahead2 - behind2
            // difference crosses zero at an inflection-point, so its sign changes
            if ((diffWas * diffNow < 0) && ((stamp[l] - then) > 100)) {
                let maybe2 = new Limit
                maybe2.value = ahead2 + behind2 //
                maybe2.time = stamp[l]
                then = stamp[l]
                yLimits.push(maybe2)
            }
        }

        diffNow = zData[2]
        diffWas = diffNow
        then = stamp[0] - 20
        for (let m = 2; m <= stamp.length - 3; m++) {
            let behind3 = zData[m] + zData[m - 1] + zData[m - 2]
            let ahead3 = zData[m] + zData[m + 1] + zData[m + 2]
            diffWas = diffNow
            diffNow = ahead3 - behind3
            // difference crosses zero at an inflection-point, so its sign changes
            if ((diffWas * diffNow < 0) && ((stamp[m] - then) > 100)) {
                let maybe3 = new Limit
                maybe3.value = ahead3 + behind3 //
                maybe3.time = stamp[m]
                zLimits.push(maybe3)
            }
        }
        // we've now finished with the raw scan data
        if (!testing) {
            stamp = []
            yData = []
            xData = []
            zData = []
        }

        // To guarantee an example of a maximum and a minimum in each dimension, we
        // need to hace scanned at least one and a bit full revolutions of the buggy...
        if ((xLimits.length < 2) || (yLimits.length < 2) || (zLimits.length < 2)) {
            return -2  // "NOT ENOUGH SCAN ROTATION"
        }

        // characterise each dim for normalisation purposes
        let [xAmp, xOff] = getRange(xLimits)
        let [yAmp, yOff] = getRange(yLimits)
        let [zAmp, zOff] = getRange(zLimits)

        // choose the two dimensions showing the biggest amplitude
        if ((xAmp > yAmp) && (xAmp > zAmp)) {  // X is the major dim
            uDim = Dimension.X
            uAmp = xAmp
            uOff2 = xOff
            uLimits = xLimits // for periodicity calculations later
            if (yAmp > zAmp) {  		      // Y is the minor dim
                vDim = Dimension.Y
                vAmp = yAmp
                vOff2 = yOff
                vLimits = yLimits
                plane = "XY"
            } else {  		                  // Z is the minor dim
                vDim = Dimension.Z
                vAmp = zAmp
                vOff2 = zOff
                vLimits = xLimits
                plane = "XZ"
            }
        }

        if ((yAmp > xAmp) && (yAmp > zAmp)) {  // Y is the major dim
            uDim = Dimension.Y
            uAmp = yAmp
            uOff2 = yOff
            uLimits = yLimits
            if (xAmp > zAmp) {  		      // X is the minor dim
                vDim = Dimension.X
                vAmp = xAmp
                vOff2 = xOff
                vLimits = xLimits
                plane = "YX"
            } else {  		                  // Z is the minor dim
                vDim = Dimension.Z
                vAmp = zAmp
                vOff2 = zOff
                vLimits = zLimits
                plane = "YZ"
            }
        }

        if ((zAmp > xAmp) && (zAmp > yAmp)) {  // Z is the major dim
            uDim = Dimension.Z
            uAmp = zAmp
            uOff2 = zOff
            uLimits = zLimits
            if (yAmp > xAmp) {  		      // Y is the minor dim
                vDim = Dimension.Y
                vAmp = yAmp
                vOff2 = yOff
                vLimits = yLimits
                plane = "ZY"
            } else {  		                  // X is the minor dim
                vDim = Dimension.X
                vAmp = xAmp
                vOff2 = xOff
                vLimits = xLimits
                plane = "ZX"
            }
        }

        // From now on, we work with the projection of the magnetic flux vector 
        // onto the selected [U,V] plane.
        // Compass-headings will be based on the strongest component (U),
        // so North equates to a positive maximum in U and a near-zero value of V.
        // Set up the normalisation factor to balance U & V projection components,
        // so that North-East equates to equal positive values of U & V
        magScale = uAmp / vAmp

        // Use the quadrant phasing of the limits[] vectors to check the sequencing 
        // of N-E-S-W compass points. If the magnetometer is mounted upside-down,
        // it will sense flux rotating in the opposite sense to the buggy's actual rotation.
        // (The scan data obvioulsy also depends on that original rotation direction)
        //          uLimits:     + 0 - 0 + 0 - 0        + 0 - 0 + 0 - 0         
        //          vLimits:     0 - 0 + 0 - 0 +        0 + 0 - 0 + 0 -
        //          heading:     N E S W N E S W        N W S E N W S E 
        // Scanned clockwise?	   matches buggy          reversed sense
        // Scanned anticlockwise?  reversed sense         matches buggy

        // see if uDim leads or trails vDim (normally we'd expect to encounter limits with opposite signs)
        if (uLimits[0].value * vLimits[0].value < 0) {
            reversed = !clockwise
        } else {
            reversed = clockwise
        }


        // Each limits vector comprises alternate hi & lo values, half a cycle apart
        // Add the averaged semi-cycles from each dimension to get an average full period
        let uSpans = uLimits.length - 1
        let vSpans = vLimits.length - 1
        let uSemi = (uLimits[uSpans].time - uLimits[0].time) / uSpans
        let vSemi = (vLimits[vSpans].time - vLimits[0].time) / vSpans
        let period = uSemi + vSemi


        // return RPM of original scan    
        return 60000 / period

    }

    // read the magnetometer and return the current heading in degrees from North
    export function degrees(): number {
        let u = 0
        let v = 0
        if (testing) {
            u = zData[test] - uOff
            v = xData[test] - vOff
            test++
            if (test == zData.length) test = 0 // roll round
        } else {
            let u = input.magneticForce(uDim) - uOff
            let v = input.magneticForce(vDim) - vOff
        }
        let val = 57.29578 * Math.atan2(u, v * magScale)
        if (reversed) {
            val = 360 - val // sensor upside-down
        }

        datalogger.log(datalogger.createCV("u", u),
            datalogger.createCV("v", v),
            datalogger.createCV("val", val))
        return (val + 360) % 360
    }

    /*export function dumpData() {
        datalogger.deleteLog()
        datalogger.setColumnTitles("t","x","y","z")
        for (let i = 0; i < stamp.length; i++) {
            datalogger.log(datalogger.createCV("t", stamp[i]),
                            datalogger.createCV("x", xData[i]),
                            datalogger.createCV("y", yData[i]),
                            datalogger.createCV("z", zData[i]))
        }
    */

    export function dumpLimits() {
        datalogger.deleteLog()
        datalogger.mirrorToSerial(true)
        datalogger.setColumnTitles("ut", "u", "vt", "v")
        for (let o = 0; o < stamp.length; o++) {
            datalogger.log(datalogger.createCV("ut", uLimits[o].time),
                datalogger.createCV("u", uLimits[o].value),
                datalogger.createCV("vt", vLimits[o].time),
                datalogger.createCV("v", vLimits[o].value))
        }
    }
}
