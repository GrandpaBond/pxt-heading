// tests go here; this will not be compiled when this package is used as an extension.
enum Config {
    Buggy,
    Jig,
    Test
}
enum Task {
    Scan,
    SetNorth,
    Measure
}


input.onButtonPressed(Button.A, function () {
    switch (nextTask) {
        case Task.Scan:
            basic.showString("S")
            basic.pause(1000)
            switch (config) {
            case Config.Buggy:
                heading.testMode(false) 
                Kitronik_Move_Motor.spin(Kitronik_Move_Motor.SpinDirections.Right, 30)
                heading.scan(4000)
                Kitronik_Move_Motor.stop()
                break
            case Config.Test:
                heading.testMode(true)
                heading.scan(1000)
                break  
            case Config.Jig:
                heading.testMode(false)
                basic.showString("?") // manually rotate jig (SMOOOOTHLY!)
                heading.scan(15000)
                music.setVolume(255)
                music.tonePlayable(2000, 500)
                basic.pause(1000)
                break
            }
            basic.showIcon(IconNames.Yes)
            basic.pause(1000)
            basic.showArrow(ArrowNames.West)
            nextTask = Task.SetNorth
            break

        case Task.SetNorth:
            basic.showString("N")
            basic.pause(1000)
            basic.clearScreen()
            spinRPM = heading.setNorth()
            basic.showNumber(Math.floor(spinRPM))
            turn30 = 60000 / (12 * spinRPM) // time needed to turn 30 degrees
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

        case Task.Measure: // Button A resets everything
            basic.showIcon(IconNames.No)
            basic.pause(1000)
            basic.clearScreen()
            basic.showArrow(ArrowNames.West)
            nextTask = Task.Scan // restart with a scan
            break
    }

})

input.onButtonPressed(Button.B, function () {
    switch (nextTask) {
        case Task.Scan: // use button A to do a scan first
            for (let i = 0; i < 5; i++) {
                basic.clearScreen()
                basic.pause(100)
                basic.showArrow(ArrowNames.West)
            }
            break

        case Task.SetNorth: // log our scan, then use button A to setNorth
            basic.showString("L")
            basic.pause(1000)
            heading.dumpData()
            basic.clearScreen()
            basic.pause(100)
            basic.showArrow(ArrowNames.West)
            break

        case Task.Measure: // OK take a new heading measurement
            if (config == Config.Buggy) {  
                basic.pause(1000)
                Kitronik_Move_Motor.spin(Kitronik_Move_Motor.SpinDirections.Right, 30)
                basic.pause(turn30) // spin to next angle
                Kitronik_Move_Motor.stop()
            } 
            // else we're testing; or we will have manually moved Jig to next test-angle...
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
            break
    }

})


input.onLogoEvent(TouchButtonEvent.Pressed, function () {
    basic.showIcon(IconNames.No)
    basic.pause(500)
    basic.clearScreen()
    switch(config) {
        case Config.Buggy:
            config = Config.Test
            basic.showString("T") // use sample data while debugging...
            heading.testMode(true) 
            break
        case Config.Test:
            config = Config.Jig
            basic.showString("J") // no buggy, but use live magnetometer
            heading.testMode(false) 
            break
        case Config.Jig:
            config = Config.Buggy
            basic.showString("B")  // normal live operation
            heading.testMode(false)
            break
    }
    basic.pause(1000)
    basic.clearScreen()
    basic.pause(200)
    nextTask = Task.Scan // new mode, so always start with a scan
    basic.showArrow(ArrowNames.West)
})


music.setVolume(255)
music.tonePlayable(2000, 500)
music.tonePlayable(1500, 500)
music.tonePlayable(1200, 500)
music.tonePlayable(1000, 1500)
let config = Config.Test
heading.testMode(true)
let nextTask = Task.Scan
basic.showArrow(ArrowNames.West)
let spinRPM = 0
let turn30 = 0
