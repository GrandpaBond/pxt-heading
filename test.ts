// tests go here; this will not be compiled when this package is used as an extension.

enum Task {
    Scan,
    SetNorth,
    Measure
}

// NOTE: check in "pxt-heading.ts" that the required test dataset is available in simulateScan()!
const dataset = "blup70_0714_1743"

function pressedA() {
    let result = 0
    switch (nextTask) {
        case Task.Scan:
            let scanTime = 6000 // ...to MANUALLY rotate turntable jig twice (SMOOOOTHLY!)
            basic.showString("S")
            basic.pause(1000)
            basic.showString("_")
            let result = heading.scanClockwise(scanTime)
            if (result == 0) {
                basic.showIcon(IconNames.Yes)
                basic.pause(1000)
                basic.clearScreen()
                basic.showArrow(ArrowNames.West)
                nextTask = Task.SetNorth
            } else {
                basic.showIcon(IconNames.Skull) // problem with scan data analysis
                basic.pause(1000)
                basic.showNumber(result)
                basic.pause(1000)
                basic.clearScreen()
                basic.showArrow(ArrowNames.West)
                nextTask = Task.Scan // try for another scan
            }
            break

        case Task.SetNorth:
            basic.showString("N")
            basic.pause(1000)
            basic.clearScreen()
            result = heading.setNorth()

            spinRPM = heading.spinRate() // ...just out of interest
            basic.showNumber(Math.floor(spinRPM))
            basic.pause(1000)
            basic.showIcon(IconNames.Yes)
            basic.pause(500)
            basic.showLeds(`
                # # . # #
                # . . . #
                . . # . .
                # . . . #
                # # . # #
                `)
            basic.pause(500)
            basic.showArrow(ArrowNames.East)
            nextTask = Task.Measure
            break

        case Task.Measure: // Button A allows new North setting
            basic.showIcon(IconNames.No)
            basic.pause(1000)
            basic.clearScreen()
            basic.showArrow(ArrowNames.West)
            nextTask = Task.SetNorth // reset new North
            break
    }
}

function pressedB() {
    switch (nextTask) {
       
        case Task.SetNorth: // sequence error!
        case Task.Scan: // use button A to do a scan first
            for (let i = 0; i < 5; i++) {
                basic.clearScreen()
                basic.pause(100)
                basic.showArrow(ArrowNames.West)
            }
            break

        case Task.Measure: // OK, take a new heading measurement
            basic.pause(200)
            basic.clearScreen()
            basic.pause(50)
            let compass = heading.degrees()
            basic.showNumber(Math.floor(compass))
            basic.pause(500)
            // now MANUALLY move to next test-angle...
            basic.showLeds(`
                    # # . # #
                    # . . . #
                    . . . . .
                    # . . . #
                    # # . # #
                    `)
            basic.pause(200)
            break
    }

}

// toggle configuration 
function toggleDebug() {
    basic.showIcon(IconNames.No)
    basic.pause(500)
    basic.clearScreen()
    if (heading.debugMode) {  // switch to Live mode
            heading.debugMode = false
            heading.scanData = []
            heading.scanTimes = []
            basic.showString("L")
    } else { // switch to Debug mode
            heading.debugMode = true
            simulateScan(dataset) // use sample data while debugging...
            basic.showString("D") 
    }

    basic.pause(1000)
    basic.clearScreen()
    basic.pause(200)
    nextTask = Task.Scan // new mode, so always start with a scan
    basic.showArrow(ArrowNames.West)
}


input.onButtonPressed(Button.A, function () {
    pressedA()
})

input.onButtonPressed(Button.B, function () {
    pressedB()
})
input.onButtonPressed(Button.AB, function () {
    toggleDebug()
})

function simulateScan(dataset: string) {
    let xData: number[] = []
    let yData: number[] = []
    let zData: number[] = []
    let xTest: number[] = []
    let yTest: number[] = []
    let zTest: number[] = []
    switch (dataset) {
        case "blup70": // angled with bottom-left upwards; dip = 70
            heading.scanTimes = [31881, 31929, 31945, 31961, 31977, 31993, 32009, 32025, 32081, 32097, 32113, 32129, 32145, 32161, 32177, 32193, 32209, 32265, 32281, 32297, 32313, 32329, 32345, 32361, 32377, 32393, 32449, 32465, 32481, 32497, 32513, 32529, 32545, 32561, 32577, 32661, 32677, 32693, 32709, 32725, 32741, 32757, 32773, 32789, 32845, 32861, 32877, 32893, 32909, 32925, 32941, 32957, 32973, 33029, 33045, 33061, 33077, 33093, 33109, 33125, 33141, 33157, 33217, 33233, 33249, 33265, 33281, 33297, 33313, 33329, 33413, 33429, 33445, 33461, 33477, 33493, 33509, 33525, 33541, 33601, 33617, 33633, 33649, 33665, 33681, 33697, 33713, 33769, 33785, 33801, 33817, 33833, 33849, 33865, 33881, 33937, 33953, 33969, 33985, 34001, 34017, 34033, 34049, 34065, 34149, 34165, 34181, 34201, 34217, 34233, 34249, 34265, 34325, 34341, 34357, 34373, 34389, 34405, 34421, 34437, 34453, 34509, 34525, 34541, 34557, 34573, 34589, 34605, 34621, 34677, 34693, 34709, 34725, 34741, 34757, 34773, 34789, 34873, 34889, 34905, 34921, 34937, 34953, 34969, 34985, 35001, 35057, 35073, 35089, 35105, 35121, 35137, 35153, 35169, 35229, 35245, 35261, 35277, 35293, 35309, 35325, 35341, 35397, 35413, 35429, 35445, 35461, 35477, 35493, 35509, 35525, 35613, 35629, 35645, 35661, 35677, 35693, 35713, 35729, 35749, 35805, 35821, 35837, 35853, 35869, 35885, 35901, 35917, 35933, 35989, 36005, 36025, 36041, 36057, 36073, 36089, 36105, 36121, 36181, 36197, 36213, 36229, 36249, 36265, 36281, 36297, 36385, 36401, 36417, 36433, 36449, 36465, 36481, 36501, 36517, 36577, 36593, 36609, 36625, 36641, 36657, 36673, 36689, 36709, 36765, 36781, 36797, 36813, 36829, 36845, 36861, 36877, 36897, 36953, 36969, 36985, 37001, 37017, 37033, 37049, 37065, 37085, 37169, 37185, 37201, 37217, 37233, 37249, 37265, 37285, 37301, 37361, 37377, 37393, 37409, 37425, 37441, 37461, 37477, 37537, 37553, 37569, 37585, 37601, 37617, 37637, 37653, 37713]
            xData = [5.61, 5.99, 6.1, 6.24, 6.38, 6.48, 6.58, 6.65, 6.86, 6.91, 6.94, 6.97, 6.99, 7.02, 7.02, 6.97, 6.94, 6.76, 6.66, 6.6, 6.53, 6.46, 6.4, 6.33, 6.25, 6.16, 5.69, 5.55, 5.45, 5.31, 5.1, 4.92, 4.75, 4.57, 4.38, 3.22, 2.98, 2.76, 2.58, 2.43, 2.25, 2.02, 1.78, 1.52, 0.57, 0.3, 0.05, -0.2, -0.49, -0.84, -1.18, -1.48, -1.76, -2.85, -3.21, -3.62, -4.01, -4.36, -4.73, -5.08, -5.4, -5.78, -7.32, -7.69, -8.02, -8.39, -8.77, -9.14, -9.54, -9.99, -12.4, -12.8, -13.17, -13.51, -13.85, -14.22, -14.6, -14.94, -15.24, -16.55, -16.93, -17.3, -17.66, -18.03, -18.37, -18.67, -18.94, -20.03, -20.35, -20.66, -20.95, -21.22, -21.48, -21.75, -22.06, -23.01, -23.23, -23.44, -23.63, -23.82, -23.99, -24.15, -24.39, -24.6, -25.37, -25.48, -25.61, -25.77, -25.83, -25.87, -25.89, -25.94, -26.1, -26.12, -26.14, -26.14, -26.13, -26.08, -26, -25.94, -25.88, -25.66, -25.62, -25.53, -25.35, -25.15, -24.95, -24.75, -24.58, -23.88, -23.68, -23.47, -23.25, -23.07, -22.84, -22.58, -22.32, -20.76, -20.4, -20.01, -19.64, -19.27, -18.96, -18.64, -18.26, -17.86, -16.45, -16, -15.54, -15.12, -14.71, -14.31, -13.88, -13.45, -11.93, -11.58, -11.21, -10.83, -10.41, -10.02, -9.64, -9.25, -7.92, -7.55, -7.18, -6.81, -6.5, -6.21, -5.92, -5.59, -5.26, -3.36, -2.98, -2.61, -2.25, -1.88, -1.52, -1.07, -0.72, -0.3, 0.78, 1.08, 1.43, 1.79, 2.08, 2.32, 2.59, 2.87, 3.1, 3.81, 4.02, 4.31, 4.56, 4.81, 5.01, 5.18, 5.35, 5.5, 5.94, 6.03, 6.13, 6.19, 6.24, 6.28, 6.31, 6.31, 6.05, 5.99, 5.92, 5.85, 5.8, 5.75, 5.68, 5.56, 5.45, 5, 4.83, 4.64, 4.42, 4.16, 3.91, 3.72, 3.49, 3.22, 2.38, 2.11, 1.85, 1.57, 1.3, 1.03, 0.76, 0.51, 0.19, -0.83, -1.16, -1.46, -1.77, -2.08, -2.41, -2.71, -2.97, -3.32, -4.98, -5.31, -5.64, -5.95, -6.26, -6.52, -6.78, -7.18, -7.51, -8.79, -9.12, -9.42, -9.77, -10.13, -10.42, -10.8, -11.15, -12.33, -12.62, -12.9, -13.18, -13.47, -13.78, -14.14, -14.41, -15.47]
            yData = [48.69, 47.65, 47.33, 47.04, 46.68, 46.24, 45.78, 45.37, 43.94, 43.53, 43.13, 42.72, 42.33, 41.95, 41.54, 41.11, 40.67, 39.09, 38.66, 38.26, 37.85, 37.45, 37.04, 36.61, 36.2, 35.83, 34.61, 34.21, 33.79, 33.4, 33.04, 32.7, 32.36, 32.05, 31.74, 30, 29.65, 29.29, 28.99, 28.74, 28.48, 28.19, 27.86, 27.53, 26.5, 26.25, 25.97, 25.65, 25.37, 25.08, 24.76, 24.49, 24.21, 23.3, 23.09, 22.89, 22.68, 22.5, 22.33, 22.12, 21.91, 21.72, 21.05, 20.91, 20.8, 20.68, 20.54, 20.37, 20.22, 20.12, 19.75, 19.71, 19.68, 19.63, 19.61, 19.6, 19.59, 19.59, 19.61, 19.84, 19.94, 20.08, 20.24, 20.38, 20.47, 20.55, 20.67, 21.17, 21.33, 21.49, 21.64, 21.81, 22, 22.15, 22.34, 23.12, 23.34, 23.56, 23.75, 24, 24.3, 24.59, 24.85, 25.14, 26.87, 27.21, 27.55, 27.98, 28.32, 28.65, 28.99, 29.32, 30.62, 30.99, 31.37, 31.74, 32.11, 32.49, 32.85, 33.22, 33.63, 35.08, 35.49, 35.91, 36.36, 36.8, 37.26, 37.75, 38.18, 39.57, 39.97, 40.36, 40.77, 41.18, 41.55, 41.92, 42.34, 44.26, 44.63, 45, 45.32, 45.67, 46.03, 46.34, 46.68, 47.04, 48.06, 48.35, 48.65, 48.94, 49.16, 49.37, 49.61, 49.83, 50.48, 50.62, 50.75, 50.9, 51.08, 51.23, 51.37, 51.5, 51.8, 51.84, 51.94, 52.03, 52.05, 52.07, 52.07, 52.05, 52.03, 51.85, 51.83, 51.8, 51.7, 51.59, 51.47, 51.31, 51.19, 51.03, 50.39, 50.15, 49.88, 49.64, 49.41, 49.18, 48.96, 48.67, 48.33, 47.07, 46.68, 46.18, 45.75, 45.33, 44.95, 44.56, 44.16, 43.77, 42.35, 41.93, 41.47, 41.06, 40.61, 40.22, 39.78, 39.36, 37.26, 36.9, 36.52, 36.14, 35.77, 35.4, 35.01, 34.54, 34.18, 32.78, 32.41, 32.02, 31.57, 31.14, 30.75, 30.39, 30.04, 29.57, 28.43, 28.15, 27.83, 27.5, 27.18, 26.85, 26.53, 26.22, 25.85, 24.82, 24.55, 24.26, 23.99, 23.78, 23.53, 23.27, 23.08, 22.85, 21.83, 21.67, 21.54, 21.38, 21.22, 21.1, 20.98, 20.81, 20.67, 20.17, 20.1, 20.07, 19.98, 19.85, 19.71, 19.59, 19.51, 19.32, 19.3, 19.28, 19.28, 19.28, 19.26, 19.22, 19.24, 19.33]
            zData = [61.37, 61.87, 62.04, 62.23, 62.42, 62.63, 62.83, 63.06, 63.86, 64.03, 64.23, 64.42, 64.55, 64.71, 64.92, 65.12, 65.28, 65.85, 66.01, 66.19, 66.37, 66.5, 66.62, 66.79, 67, 67.17, 67.57, 67.6, 67.65, 67.76, 67.85, 67.98, 68.09, 68.18, 68.28, 68.64, 68.7, 68.69, 68.7, 68.76, 68.85, 68.93, 68.96, 69, 69.08, 69.06, 69.04, 69.07, 69.12, 69.1, 69.04, 68.99, 68.97, 68.93, 68.88, 68.83, 68.78, 68.72, 68.67, 68.59, 68.53, 68.45, 68.14, 68.05, 67.94, 67.83, 67.72, 67.64, 67.51, 67.37, 66.72, 66.56, 66.36, 66.15, 65.99, 65.86, 65.71, 65.54, 65.35, 64.65, 64.45, 64.26, 64.04, 63.83, 63.63, 63.42, 63.21, 62.46, 62.23, 62, 61.79, 61.57, 61.36, 61.17, 60.99, 60.25, 60.03, 59.82, 59.6, 59.4, 59.19, 58.95, 58.76, 58.58, 57.42, 57.2, 56.95, 56.69, 56.55, 56.37, 56.17, 56.02, 55.48, 55.28, 55.09, 54.95, 54.78, 54.57, 54.38, 54.24, 54.13, 53.76, 53.6, 53.44, 53.29, 53.14, 53.02, 52.87, 52.73, 52.28, 52.17, 52.06, 51.99, 51.96, 51.92, 51.87, 51.8, 51.47, 51.46, 51.49, 51.5, 51.49, 51.49, 51.52, 51.59, 51.67, 51.91, 51.98, 52.04, 52.07, 52.09, 52.15, 52.23, 52.3, 52.54, 52.61, 52.67, 52.75, 52.9, 53.11, 53.31, 53.45, 53.9, 54.02, 54.11, 54.21, 54.39, 54.59, 54.74, 54.9, 55.08, 56.06, 56.23, 56.4, 56.62, 56.83, 57.03, 57.32, 57.57, 57.79, 58.33, 58.53, 58.74, 58.97, 59.24, 59.53, 59.77, 59.95, 60.14, 61.02, 61.32, 61.66, 61.92, 62.22, 62.5, 62.76, 63, 63.26, 64.16, 64.32, 64.43, 64.61, 64.92, 65.13, 65.31, 65.49, 66.29, 66.43, 66.59, 66.76, 66.91, 67.05, 67.15, 67.26, 67.37, 67.75, 67.89, 68.03, 68.14, 68.26, 68.35, 68.4, 68.44, 68.51, 68.71, 68.75, 68.75, 68.77, 68.82, 68.87, 68.94, 69.02, 69.07, 69.11, 69.1, 69.08, 69.05, 69.02, 68.98, 68.98, 69, 68.97, 68.57, 68.5, 68.46, 68.38, 68.29, 68.23, 68.18, 68.13, 68.05, 67.69, 67.58, 67.45, 67.35, 67.25, 67.16, 67.06, 66.99, 66.64, 66.54, 66.44, 66.35, 66.29, 66.23, 66.1, 65.99, 65.53]
            xTest = [-25.08, -26.12, -17.66, -4.49, 4.99, 6.55, -2.46, -15.36, -25.47, -26.06, -18.31, -4.9, 5.21, 6.13, -2.34, -16.07, -24.82, -26.35, -18.15, -5.31, 4.85, 6.33, -1.89, -15.62, -25.71, -26.65, -18.33, -5, 4.72, 5.88, -2.5, -16.29, -25.92, -26.85, -18.47, -5.35, 4.68, 5.7, -2.68, -16.25, -25.69, -26.87, -18.45, -5.47, 4.68, 5.94, -2.52, -15.97, -25.31, -26.35, -18.49, -5.12, 5.05, 5.68, -2.68, -16.13, -25.8, -27.42, -18.41, -5.49, 4.76, 5.41, -2.83, -16.38, -26.04]
            yTest = [24.46, 36.27, 48.55, 52.71, 46.86, 34.78, 22.73, 18.1, 23.43, 35.59, 47.44, 52.37, 46.57, 34.56, 22.83, 18.65, 24.02, 35.91, 47.7, 52.06, 46.71, 33.93, 22.4, 17.92, 22.93, 35.22, 47.14, 51.82, 46.07, 33.57, 21.92, 17.31, 22.83, 34.17, 46.61, 51.12, 45.48, 33.55, 21.03, 17.11, 22.71, 34.46, 46.43, 51.28, 45.5, 33.17, 21.49, 16.97, 22.3, 33.79, 46.49, 51.32, 45.58, 32.84, 21.23, 16.45, 22.65, 34.32, 46.11, 50.39, 45.02, 32.68, 21.43, 16.73, 21.88]
            zTest = [58.77, 52.87, 51.14, 54.79, 61.46, 67.49, 69.05, 65.55, 58.29, 52.56, 50.62, 54.1, 61.38, 67.57, 69.12, 65.63, 58.54, 52.35, 50.72, 53.85, 60.86, 67.34, 68.74, 65.4, 58.67, 52.74, 51.22, 54.64, 61.25, 67.74, 69.49, 65.61, 58.98, 52.81, 50.95, 54.89, 61.4, 67.24, 69.39, 65.57, 58.54, 52.45, 51.43, 54.79, 61.34, 67.8, 69.47, 65.88, 58.4, 52.93, 51.7, 55.41, 61.38, 67.34, 69.8, 65.8, 58.75, 52.77, 50.99, 54.6, 61.69, 67.36, 69.39, 65.68, 58.42]
            break


        case "blup70_0714_1743":
            heading.scanTimes = [32009, 32057, 32073, 32089, 32105, 32121, 32137, 32193, 32209, 32225, 32241, 32257, 32273, 32289, 32305, 32361, 32377, 32393, 32409, 32425, 32441, 32457, 32473, 32529, 32545, 32561, 32577, 32593, 32609, 32625, 32713, 32729, 32745, 32761, 32777, 32793, 32809, 32825, 32885, 32901, 32917, 32933, 32949, 32965, 32981, 33037, 33053, 33069, 33085, 33101, 33117, 33133, 33149, 33205, 33221, 33237, 33253, 33269, 33285, 33301, 33385, 33401, 33417, 33433, 33449, 33465, 33481, 33497, 33553, 33569, 33585, 33601, 33617, 33633, 33649, 33665, 33721, 33737, 33753, 33769, 33785, 33801, 33817, 33873, 33889, 33905, 33921, 33937, 33953, 33969, 33985, 34069, 34085, 34101, 34117, 34133, 34149, 34165, 34193, 34225, 34241, 34257, 34273, 34289, 34305, 34321, 34381, 34397, 34413, 34429, 34445, 34461, 34477, 34493, 34549, 34565, 34581, 34597, 34613, 34629, 34645, 34729, 34745, 34761, 34777, 34793, 34809, 34825, 34841, 34897, 34913, 34929, 34945, 34961, 34977, 34993, 35049, 35065, 35081, 35097, 35113, 35129, 35145, 35161, 35217, 35233, 35249, 35265, 35281, 35297, 35313, 35329, 35413, 35429, 35445, 35461, 35477, 35493, 35509, 35565, 35581, 35597, 35613, 35629, 35645, 35661, 35677, 35733, 35749, 35765, 35781, 35797, 35813, 35829, 35889, 35905, 35921, 35941, 35957, 35977, 35993, 36009, 36093, 36109, 36125, 36141, 36157, 36173, 36189, 36205, 36265, 36285, 36301, 36317, 36333, 36349, 36365, 36425, 36441, 36457, 36473, 36489, 36509, 36525, 36541, 36601, 36617, 36633, 36649, 36665, 36681, 36697, 36717, 36801, 36817, 36833, 36849, 36865, 36881, 36897, 36957, 36973, 36993, 37009, 37025, 37041, 37057, 37073, 37133, 37149, 37165, 37185, 37201, 37217, 37233, 37293, 37309, 37325, 37341, 37357, 37377, 37393, 37409, 37513, 37529, 37545, 37561, 37577, 37597, 37613, 37629, 37689, 37705, 37721, 37737, 37753, 37773, 37789, 37849]
            xData = [887.59, 889.13, 889.71, 890.29, 890.92, 891.59, 892.27, 894.78, 895.5, 896.22, 896.95, 897.69, 898.51, 899.29, 899.92, 902.1, 902.74, 903.37, 903.98, 904.58, 905.16, 905.7, 906.19, 907.75, 908.14, 908.47, 908.75, 909.03, 909.29, 909.5, 910.21, 910.23, 910.19, 910.13, 910.01, 909.81, 909.59, 909.38, 908.27, 907.95, 907.63, 907.23, 906.74, 906.24, 905.79, 904.17, 903.63, 903.04, 902.39, 901.71, 901.1, 900.47, 899.78, 897.25, 896.48, 895.71, 894.94, 894.18, 893.42, 892.63, 888.97, 888.33, 887.69, 887.06, 886.45, 885.88, 885.39, 884.91, 883.48, 883.16, 882.84, 882.57, 882.36, 882.23, 882.15, 882.08, 881.96, 881.99, 882.13, 882.31, 882.48, 882.68, 882.87, 883.82, 884.23, 884.68, 885.14, 885.63, 886.13, 886.65, 887.19, 890.34, 891.01, 891.75, 892.55, 893.35, 894.14, 894.94, 896.39, 898.01, 898.79, 899.56, 900.32, 901.05, 901.8, 902.54, 904.97, 905.53, 906.04, 906.5, 906.89, 907.31, 907.76, 908.13, 909.11, 909.33, 909.49, 909.6, 909.69, 909.74, 909.72, 909.25, 909.13, 908.98, 908.75, 908.48, 908.18, 907.89, 907.58, 906.27, 905.85, 905.39, 904.88, 904.39, 903.86, 903.32, 901.4, 900.8, 900.13, 899.43, 898.78, 898.14, 897.5, 896.82, 894.41, 893.76, 893.12, 892.48, 891.85, 891.17, 890.49, 889.84, 886.85, 886.33, 885.85, 885.4, 884.95, 884.48, 884.01, 882.74, 882.48, 882.23, 882.03, 881.87, 881.7, 881.53, 881.44, 881.36, 881.36, 881.42, 881.55, 881.72, 881.89, 882.06, 883.09, 883.45, 883.83, 884.43, 884.98, 885.72, 886.35, 886.94, 890.46, 891.23, 891.97, 892.68, 893.45, 894.29, 895.16, 896.04, 899.26, 900.29, 901.1, 901.88, 902.62, 903.3, 903.95, 906.16, 906.66, 907.11, 907.53, 907.91, 908.32, 908.62, 908.92, 909.63, 909.73, 909.74, 909.7, 909.68, 909.67, 909.62, 909.46, 908.19, 907.91, 907.64, 907.3, 906.94, 906.64, 906.29, 904.66, 904.21, 903.61, 903.11, 902.59, 902.04, 901.47, 900.89, 898.83, 898.25, 897.66, 896.93, 896.36, 895.78, 895.17, 892.83, 892.24, 891.69, 891.13, 890.52, 889.72, 889.11, 888.53, 885.11, 884.67, 884.28, 883.88, 883.49, 883.08, 882.74, 882.42, 881.65, 881.51, 881.38, 881.31, 881.27, 881.28, 881.32, 881.75, 881.16]
            yData = [1586.86, 1587.98, 1588.36, 1588.68, 1589, 1589.33, 1589.65, 1590.42, 1590.58, 1590.69, 1590.77, 1590.85, 1590.91, 1590.88, 1590.8, 1590.45, 1590.3, 1590.09, 1589.86, 1589.56, 1589.22, 1588.9, 1588.58, 1587.2, 1586.76, 1586.24, 1585.67, 1585.17, 1584.69, 1584.12, 1580.57, 1579.91, 1579.24, 1578.54, 1577.81, 1577.08, 1576.37, 1575.7, 1573.32, 1572.64, 1571.95, 1571.27, 1570.54, 1569.85, 1569.19, 1567.16, 1566.65, 1566.15, 1565.65, 1565.18, 1564.74, 1564.34, 1563.96, 1562.9, 1562.71, 1562.58, 1562.46, 1562.37, 1562.32, 1562.28, 1562.83, 1563.04, 1563.27, 1563.55, 1563.85, 1564.18, 1564.57, 1564.97, 1566.67, 1567.26, 1567.86, 1568.43, 1569.04, 1569.67, 1570.34, 1571.08, 1573.56, 1574.22, 1574.9, 1575.69, 1576.44, 1577.15, 1577.85, 1580.1, 1580.72, 1581.35, 1581.96, 1582.56, 1583.19, 1583.75, 1584.3, 1586.89, 1587.3, 1587.66, 1587.97, 1588.28, 1588.58, 1588.85, 1589.2, 1589.44, 1589.49, 1589.44, 1589.33, 1589.18, 1588.97, 1588.75, 1587.81, 1587.52, 1587.15, 1586.71, 1586.19, 1585.66, 1585.17, 1584.7, 1582.82, 1582.22, 1581.61, 1581.01, 1580.42, 1579.8, 1579.15, 1575.98, 1575.36, 1574.75, 1574.14, 1573.54, 1572.96, 1572.38, 1571.81, 1569.68, 1569.09, 1568.52, 1567.99, 1567.52, 1567.02, 1566.5, 1564.87, 1564.47, 1564.16, 1563.89, 1563.54, 1563.15, 1562.83, 1562.58, 1562.06, 1561.96, 1561.89, 1561.84, 1561.78, 1561.75, 1561.81, 1561.93, 1562.98, 1563.29, 1563.62, 1563.94, 1564.28, 1564.64, 1565.02, 1566.62, 1567.08, 1567.55, 1568.05, 1568.58, 1569.14, 1569.68, 1570.18, 1572.16, 1572.79, 1573.48, 1574.17, 1574.85, 1575.52, 1576.18, 1578.69, 1579.38, 1580.11, 1581.01, 1581.7, 1582.52, 1583.15, 1583.72, 1586.25, 1586.71, 1587.18, 1587.59, 1587.97, 1588.29, 1588.54, 1588.75, 1589.03, 1588.98, 1588.89, 1588.76, 1588.61, 1588.43, 1588.2, 1586.94, 1586.51, 1586.08, 1585.67, 1585.27, 1584.71, 1584.17, 1583.53, 1581.21, 1580.6, 1579.94, 1579.27, 1578.61, 1577.96, 1577.31, 1576.5, 1573.1, 1572.49, 1571.93, 1571.34, 1570.65, 1569.96, 1569.34, 1567.43, 1566.95, 1566.37, 1565.97, 1565.57, 1565.15, 1564.73, 1564.34, 1563.05, 1562.75, 1562.5, 1562.25, 1562.04, 1561.85, 1561.7, 1561.32, 1561.29, 1561.28, 1561.3, 1561.35, 1561.48, 1561.58, 1561.71, 1563.58, 1564, 1564.44, 1564.89, 1565.37, 1565.9, 1566.32, 1566.79, 1568.82, 1569.37, 1569.94, 1570.58, 1571.21, 1571.97, 1572.6, 1575.1, 1566.09]
            zData = [424.65, 424.91, 425.05, 425.15, 425.24, 425.37, 425.57, 426.47, 426.72, 426.93, 427.14, 427.37, 427.62, 427.9, 428.19, 429.33, 429.66, 429.95, 430.25, 430.57, 430.88, 431.19, 431.54, 432.76, 433.13, 433.54, 433.92, 434.26, 434.56, 434.86, 436.42, 436.7, 436.97, 437.25, 437.52, 437.77, 437.99, 438.19, 438.85, 439.02, 439.16, 439.28, 439.38, 439.46, 439.57, 439.72, 439.7, 439.67, 439.6, 439.47, 439.38, 439.35, 439.31, 438.83, 438.62, 438.39, 438.15, 437.87, 437.58, 437.27, 435.39, 435.02, 434.7, 434.36, 434, 433.62, 433.21, 432.78, 431.34, 430.96, 430.6, 430.24, 429.84, 429.48, 429.16, 428.84, 427.79, 427.51, 427.26, 427, 426.75, 426.57, 426.44, 426.01, 425.91, 425.82, 425.7, 425.61, 425.55, 425.49, 425.45, 425.69, 425.8, 425.95, 426.12, 426.33, 426.55, 426.74, 427.16, 427.75, 428.05, 428.39, 428.81, 429.2, 429.55, 429.91, 431.31, 431.67, 432.02, 432.38, 432.72, 433.09, 433.45, 433.81, 435.05, 435.36, 435.68, 436.02, 436.32, 436.61, 436.92, 438.33, 438.53, 438.77, 438.97, 439.11, 439.29, 439.48, 439.61, 439.84, 439.93, 439.97, 439.96, 439.96, 439.95, 439.97, 439.8, 439.67, 439.55, 439.47, 439.33, 439.2, 439.08, 438.92, 438.25, 438.04, 437.84, 437.57, 437.28, 437, 436.68, 436.35, 434.66, 434.3, 433.91, 433.56, 433.23, 432.91, 432.58, 431.48, 431.21, 430.9, 430.59, 430.26, 429.91, 429.61, 429.34, 428.46, 428.2, 427.98, 427.74, 427.5, 427.33, 427.16, 426.47, 426.32, 426.18, 426, 425.91, 425.8, 425.74, 425.75, 426.1, 426.22, 426.34, 426.47, 426.66, 426.88, 427.1, 427.35, 428.54, 428.98, 429.35, 429.7, 430.04, 430.41, 430.79, 432.22, 432.61, 432.99, 433.35, 433.71, 434.16, 434.5, 434.81, 436.04, 436.39, 436.7, 436.95, 437.26, 437.6, 437.9, 438.22, 439.05, 439.15, 439.24, 439.37, 439.52, 439.68, 439.78, 439.89, 439.95, 440.02, 440.02, 439.95, 439.88, 439.85, 439.85, 439.58, 439.43, 439.29, 439.12, 438.97, 438.8, 438.59, 437.8, 437.55, 437.27, 437.06, 436.82, 436.41, 436.1, 435.79, 433.47, 433.13, 432.81, 432.47, 432.13, 431.72, 431.38, 431.05, 429.9, 429.59, 429.32, 429.04, 428.72, 428.38, 428.13, 427.27, 430.74]
            xTest = [881.04, 880.44, 889.41, 901.18, 910.09, 911.06, 901.67, 889.44, 880.74, 880.39, 888.66, 900.99, 910.05, 910.09, 901.37, 889.26, 880.29, 879.88, 888.69, 900.51, 909.99, 909.77, 901.22, 888.58, 879.79]
            yTest = [1566.06, 1577.64, 1588.18, 1591.76, 1585.86, 1573.95, 1562.89, 1559.31, 1565.21, 1576.89, 1587.86, 1591.16, 1584.79, 1573.18, 1562.72, 1559.14, 1565.25, 1576.5, 1587.41, 1590.28, 1584.86, 1572.92, 1562.46, 1558.41, 1564.11]
            zTest = [430.54, 425.51, 424.33, 428.21, 434.21, 439.48, 440.53, 437.04, 430.59, 425.85, 424.91, 428.08, 434.46, 439.11, 440.31, 436.91, 430.22, 425.46, 424.44, 427.63, 434.36, 439.48, 440.25, 436.84, 430.5]
            break
    }

    // transpose the three arrays into array of triples
    for (let i = 0; i < heading.scanTimes.length; i++) {
        let xyz = []
        xyz.push(xData[i])
        xyz.push(yData[i])
        xyz.push(zData[i])
        heading.scanData.push(xyz)
    }
    // do the same for the test cases
    for (let n = 0; n < xTest.length; n++) {
        let xyz = []
        xyz.push(xTest[n])
        xyz.push(yTest[n])
        xyz.push(zTest[n])
        heading.testData.push(xyz)
    }
}

let nextTask: Task
heading.debugMode = true // --> normal live action when A+B pressed

for (let i = 0; i < 5; i++) {
    basic.clearScreen()
    basic.pause(100)
    // invite A+B press
    basic.showLeds(`
                    . . . . .
                    . # . # .
                    # # . # #
                    . # . # .
                    . . . . .
                    `)
}
let spinRPM = 0
