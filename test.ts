// tests go here; this will not be compiled when this package is used as an extension.
enum Config {
    Live, // Normal usage (but use turntable Jig to pretend it's on a buggy)
    Capture, // Acquire new test datasets, using turntable Jig
    Debug, // Test & debug (dataset selection is preset in code below)
}
enum Task {
    Scan,
    SetNorth,
    Measure
}
// NOTe: check in pxt-heading.ts that the required test dataset is not commented-out in simulateScan()!
//const dataset = "angled"
//const dataset = "yup70"
//const dataset = "zdown70"
//const dataset = "zdown15"
//const dataset = "tldown70"
//const dataset = "tldown0"
const dataset = "dashboard70"


function performSetup() {
    let result = 0
    switch (nextTask) {
        case Task.Scan:
            let scanTime = 6000 // ...to manually rotate turntable jig twice (SMOOOOTHLY!)
            basic.showString("S")
            basic.pause(1000) 
            switch (config) {
                case Config.Live:
                    heading.setMode(Mode.Normal,"")
                    break
                case Config.Debug:
                    heading.setMode(Mode.Debug, dataset)
                    scanTime = 500 // don't wait around
                    break
                case Config.Capture:
                    heading.setMode(Mode.Capture, "")
                    break
            }
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
                basic.showIcon(IconNames.Skull) // scan failed 
                basic.pause(1000)
                basic.showNumber(result)
                basic.pause(1000)
                basic.clearScreen()
                basic.showArrow(ArrowNames.West)
                nextTask = Task.Scan // restart with a fresh scan
            } else {
                spinRPM = heading.spinRPM() // ...out of interest
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

        case Task.Measure: // Button A resets everything
            basic.showIcon(IconNames.No)
            basic.pause(1000)
            basic.clearScreen()
            basic.showArrow(ArrowNames.West)
            nextTask = Task.Scan // restart with a scan
            datalogger.deleteLog()
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
            // now manually move to next test-angle...
            basic.showLeds(`
                    # # . # #
                    # . . . #
                    . . # . .
                    # . . . #
                    # # . # #
                    `)
            basic.pause(200)
            break
    }

}

// rotate three-state configuration 
function nextConfig() {
    basic.showIcon(IconNames.No)
    basic.pause(500)
    basic.clearScreen()
    switch(config) {
        case Config.Live:
            config = Config.Debug
            heading.setMode(Mode.Debug, dataset)
            basic.showString("D") // use sample data while debugging...
            break
        case Config.Debug:
            config = Config.Capture
            basic.showString("C") // no buggy, but use live magnetometer
            heading.setMode(Mode.Capture, "")
            break
        case Config.Capture:
            config = Config.Live
            basic.showString("L")  // normal live operation
            heading.setMode(Mode.Normal, "")
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
let config = Config.Live
nextConfig() // always start with Config.Debug ...until it all works!
let spinRPM = 0
