enum Config {
    Buggy,
    Jig,
    Test
}
enum Task {
    Scan,
    Analyse,
    Measure
}

// tests go here; this will not be compiled when this package is used as an extension.
input.onButtonPressed(Button.A, function () {
    switch (task) {
        case Task.Scan:
            basic.showString("S")
            if (config == Config.Buggy) {
                heading.testMode(false) 
                Kitronik_Move_Motor.spin(Kitronik_Move_Motor.SpinDirections.Right, 30)
                heading.scan(4000)
                Kitronik_Move_Motor.stop()
            } else { // using jig, or in test mode
                heading.scan(2000)   
            }
            basic.pause(1000)
            basic.showIcon(IconNames.Yes)
            basic.pause(500)
            basic.showArrow(ArrowNames.West)
            task = Task.Analyse
        break

        case Task.Analyse:
            basic.showString("A")
            basic.pause(1000)
            basic.clearScreen()
            spinRPM = heading.analyseScan()
            basic.showNumber(spinRPM)
            turn30 = 60000 / (12 * spinRPM) // time to turn 30 degree
            datalogger.setColumnTitles("uRaw", "vRaw", "u", "v", "val")
            basic.pause(1000)
            basic.showIcon(IconNames.Yes)
            basic.pause(500)
            task = Task.Measure
            basic.showArrow(ArrowNames.East)
        break

        case Task.Measure:
            basic.showIcon(IconNames.No)
            basic.pause(1000)
            basic.clearScreen()
            task = Task.Scan // restart with a scan
            basic.showArrow(ArrowNames.West)
        break
    }

})

input.onButtonPressed(Button.B, function () {
    if (task != Task.Measure) {
        basic.clearScreen()
        basic.pause(200)
        basic.showArrow(ArrowNames.West)
        basic.pause(200)
        basic.clearScreen()
        basic.pause(200)
        basic.showArrow(ArrowNames.West)
        basic.pause(200)
        basic.clearScreen()
        basic.pause(200)
        basic.showArrow(ArrowNames.West)

    } else {
        if (config == Config.Buggy){
            Kitronik_Move_Motor.spin(Kitronik_Move_Motor.SpinDirections.Right, 30)
            basic.pause(turn30)
            Kitronik_Move_Motor.stop()
        }
        let compass = heading.degrees()
        basic.clearScreen()   
        basic.pause(500)
        basic.showNumber(Math.floor(compass))
        basic.clearScreen()   
        basic.pause(500)
        basic.showArrow(ArrowNames.East)
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
            basic.showString("J") // no buggy, but use magnetometer
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
    basic.pause(1000)
    task = Task.Scan
    basic.showArrow(ArrowNames.West)
})

let config = Config.Buggy
let task = Task.Scan
basic.showArrow(ArrowNames.West) // normal live operation
let spinRPM = 0
let turn30 = 0
