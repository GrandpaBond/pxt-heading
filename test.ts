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


input.onButtonPressed(Button.A, function () {
    performSetup()
})

function performSetup() {
    let result = 0
    //checkLogging()
    switch (nextTask) {
        case Task.Scan:
            basic.showString("S")
            basic.pause(1000)
            switch (config) {
                case Config.Buggy:
                    heading.testDataset("NONE")
                    Kitronik_Move_Motor.spin(Kitronik_Move_Motor.SpinDirections.Right, 30)
                    heading.scan(4000)
                    Kitronik_Move_Motor.stop()
                    break
                case Config.Debug:            
                    // heading.testDataset("Zdn70") // X-Axis upwards; round X-Axis; 70-degreee dip
                    // heading.testDataset("Yup70") // Y-Axis upwards; round Y-Axis; 70-degreee dip
                    // heading.testDataset("Zdn70") // Z-Axis downwards; round Z-Axis; 70-degreee dip
                    heading.testDataset("strange") // No axis aligned with vertical rotation;  70-degreee dip
                    heading.scan(1000)
                    break
                case Config.Jig:
                    heading.testDataset("NONE")
                    basic.showString("?") // manually rotate jig (SMOOOOTHLY!)
                    heading.scan(5000)
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

input.onButtonPressed(Button.AB, function () {
    nextConfig()
})

// rotate three-state configuration 
function nextConfig() {
    basic.showIcon(IconNames.No)
    basic.pause(500)
    basic.clearScreen()
    switch(config) {
        case Config.Buggy:
            config = Config.Debug
            basic.showString("D") // use sample data while debugging...
            break
        case Config.Debug:
            config = Config.Jig
            basic.showString("J") // no buggy, but use live magnetometer
            heading.testDataset("NONE")
            break
        case Config.Jig:
            config = Config.Buggy
            basic.showString("B")  // normal live operation
            heading.testDataset("NONE")
            break
    }
    basic.pause(1000)
    basic.clearScreen()
    basic.pause(200)
    nextTask = Task.Scan // new mode, so always start with a scan
    basic.showArrow(ArrowNames.West)
}


heading.setLogMode(true)
let nextTask: Task
let config = Config.Buggy
nextConfig() // start with Config.Debug ...until it all works!
let spinRPM = 0
let turn45 = 0
