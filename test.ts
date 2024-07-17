// tests go here; this will not be compiled when this package is used as an extension.
enum Task {
    Scan,
    SetNorth,
    Measure
}

function performSetup() {
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
            } else { 
                basic.showIcon(IconNames.Skull) // problem with scan data analysis
                basic.pause(1000)
                basic.showNumber(result)
                basic.pause(1000)
                basic.clearScreen()
                basic.showArrow(ArrowNames.West)
                nextTask = Task.Scan // restart with a fresh scan
            } 
            basic.clearScreen()
            basic.showArrow(ArrowNames.West)
            nextTask = Task.SetNorth
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


input.onButtonPressed(Button.A, function () {
    performSetup()
})

input.onButtonPressed(Button.B, function () {
    measure()
})

let nextTask: Task = Task.Scan
let spinRPM = 0
for(let i = 0; i < 5; i++) {
    basic.showArrow(ArrowNames.West)
    basic.pause(100)
    basic.clearScreen()
    basic.pause(100)
}
