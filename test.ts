// tests go here; this will not be compiled when this package is used as an extension.
enum Config {
    Live, // Normal usage (but use turntable Jig to pretend it's on a buggy)
    Capture, // Acquire new test datasets, using turntable Jig
    Analyse, // Test & debug (dataset selection is preset in code below)
    Trace, // Gather full diagnostics using dataLogger
}
enum Task {
    Scan,
    SetNorth,
    Measure
}
// NOTE: check in "pxt-heading.ts" that the required test dataset is available in simulateScan()!
const dataset = "blup70_0714_1743"

function performSetup() {
    let result = 0
    switch (nextTask) {
        case Task.Scan:
            let scanTime = 6000 // ...to MANUALLY rotate turntable jig twice (SMOOOOTHLY!)
            basic.showString("S")
            basic.pause(1000)
            basic.showString("_")
            heading.scanClockwise(scanTime)
            basic.showIcon(IconNames.Yes)
            basic.pause(1000)
            basic.clearScreen()
            basic.showArrow(ArrowNames.West)
            nextTask = Task.SetNorth
            break

        case Task.SetNorth:
            basic.showString("N")
            basic.pause(1000)
            basic.clearScreen()
            result = heading.setNorth()

            if (result < 0) {
                basic.showIcon(IconNames.Skull) // problem with scan data analysis
                basic.pause(1000)
                basic.showNumber(result)
                basic.pause(1000)
                basic.clearScreen()
                basic.showArrow(ArrowNames.West)
                nextTask = Task.Scan // restart with a fresh scan
            } else {
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
            }
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

function measure() {
    switch (nextTask) {
        // ? sequence error?
        case Task.SetNorth: 
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

// rotate four-state configuration 
function nextConfig() {
    basic.showIcon(IconNames.No)
    basic.pause(500)
    basic.clearScreen()
    switch(config) {
        case Config.Live:
            config = Config.Analyse
            heading.resetMode(Mode.Analyse, dataset)
            basic.showString("A") // use sample data while debugging...
            break
        case Config.Analyse:
            config = Config.Capture
            heading.resetMode(Mode.Capture, "")
            basic.showString("C") // no buggy, but use live magnetometer
            break
        case Config.Capture:
            config = Config.Trace
            heading.resetMode(Mode.Trace, dataset)
            basic.showString("T")  // gather full diagnostics
            break
        case Config.Trace:
            config = Config.Live
            heading.resetMode(Mode.Normal, "")
            basic.showString("L")  // normal live operation
            break
    }
    basic.pause(1000)
    basic.clearScreen()
    basic.pause(200)
    nextTask = Task.Scan // new mode, so always start with a scan
    basic.showArrow(ArrowNames.West)
}


input.onButtonPressed(Button.A, function () {
    performSetup()
})

input.onButtonPressed(Button.B, function () {
    measure()
})
input.onButtonPressed(Button.AB, function () {
    nextConfig()
})

let nextTask: Task
let config = Config.Capture // --> Config.Live when A+B pressed

for(let i = 0; i < 5; i++) {
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
