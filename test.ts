function performSetup() {
    switch (nextTask) {
        case Task.Scan:
            let scanTime = 5000 // ...to MANUALLY rotate twice (SMOOOOTHLY!)
            basic.showString("S")
            basic.pause(1000)
            basic.showString("_")
            result = heading.scanClockwise(scanTime)
            if (result == 0) {
                basic.showIcon(IconNames.Yes)
                nextTask = Task.SetNorth
            } else {
                basic.showIcon(IconNames.Skull) // problem with scan data analysis
                basic.pause(1000)
                basic.showNumber(result)
            }
            basic.pause(1000)
            basic.clearScreen()
            basic.showArrow(ArrowNames.West)
            break

        case Task.SetNorth:
            basic.showString("N")
            basic.pause(1000)
            basic.clearScreen()
            result = heading.setNorth()
            if (result == 0) {
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
input.onButtonPressed(Button.A, function () {
    performSetup()
})
input.onButtonPressed(Button.B, function () {
    measure()
})
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
let result = 0
let spinRPM = 0
enum Task {
    Scan,
    SetNorth,
    Measure
}
let nextTask: Task = Task.Scan
basic.showArrow(ArrowNames.West)
