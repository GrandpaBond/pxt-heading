// tests go here; this will not be compiled when this package is used as an extension.
enum Config {
    Buggy, // Normal usage, in a Kitronik Move Motor buggy
    Jig, // Acquiring new test datasets, using a special rotating jig
    Debug, // Test & debug (dataset selection is preset in code below)
}
enum Task {
    Scan,
    SetNorth,
    Measure
}

const dataset = "tldown70"


function performSetup() {
    let result = 0
    switch (nextTask) {
        case Task.Scan:
            basic.showString("S")
            basic.pause(1000)
            switch (config) {
                case Config.Buggy:
                    heading.setMode(Mode.Normal,"")
                    Kitronik_Move_Motor.spin(Kitronik_Move_Motor.SpinDirections.Right, 30)
                    heading.scan(6000)
                    Kitronik_Move_Motor.stop()
                    break
                case Config.Debug:
                    heading.setMode(Mode.Debug, dataset)
                    heading.scan(1000)
                    break
                case Config.Jig:
                    heading.setMode(Mode.Capture, "")
                    basic.showString("?") // manually rotate jig (SMOOOOTHLY!)
                    heading.scan(7000)
                    basic.pause(1000)
                    break
            }
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
                spinRPM = heading.spinRPM()
                basic.showNumber(Math.floor(spinRPM))
                turn45 = 60000 / (8 * spinRPM) // time needed to turn 45 degrees
                basic.showIcon(IconNames.Yes)
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
            basic.pause(1000)
            let compass = heading.degrees()
            basic.clearScreen()   
            basic.pause(500)
            basic.showNumber(Math.floor(compass))
            basic.clearScreen()   
            basic.pause(200)
            basic.showLeds(`
                    # # . # #
                    # . . . #
                    . . # . .
                    # . . . #
                    # # . # #
                    `)
            basic.pause(500)
            basic.showArrow(ArrowNames.East)
            // On the live buggy, move 45 degrees to next test angle;
            if (config == Config.Buggy) {  
                basic.pause(1000)
                Kitronik_Move_Motor.spin(Kitronik_Move_Motor.SpinDirections.Right, 30)
                basic.pause(turn45) // spin to next angle
                Kitronik_Move_Motor.stop()
            } 
            // else we're simulating; or must manually move Jig to next test-angle...
            break
    }

}

// rotate three-state configuration 
function nextConfig() {
    basic.showIcon(IconNames.No)
    basic.pause(500)
    basic.clearScreen()
    switch(config) {
        case Config.Buggy:
            config = Config.Debug
            heading.setMode(Mode.Debug, dataset)
            basic.showString("D") // use sample data while debugging...
            break
        case Config.Debug:
            config = Config.Jig
            basic.showString("J") // no buggy, but use live magnetometer
            heading.setMode(Mode.Capture, "")
            break
        case Config.Jig:
            config = Config.Buggy
            basic.showString("B")  // normal live operation
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
let config = Config.Buggy
nextConfig() // always start with Config.Debug ...until it all works!
let spinRPM = 0
let turn45 = 0
