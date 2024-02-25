// tests go here; this will not be compiled when this package is used as an extension.
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
                basic.showString("?")
                heading.scan(10000)
                music.setVolume(255)
                music.tonePlayable(Note.C, music.beat(BeatFraction.Whole))
                basic.pause(1000)
                break
            }
            basic.showIcon(IconNames.Yes)
            basic.pause(1000)
            basic.showArrow(ArrowNames.West)
            nextTask = Task.Analyse
        break

        case Task.Analyse:
            basic.showString("A")
            basic.pause(1000)
            basic.clearScreen()
            spinRPM = heading.analyseScan()
            basic.showNumber(Math.floor(spinRPM))
            turn30 = 60000 / (12 * spinRPM) // time to turn 30 degree
            datalogger.setColumnTitles("uRaw", "vRaw", "u", "v", "val")
            basic.pause(1000)
            basic.showIcon(IconNames.Yes)
            basic.pause(500)
            nextTask = Task.Measure
            basic.showArrow(ArrowNames.East)
        break

        case Task.Measure:
            basic.showIcon(IconNames.No)
            basic.pause(1000)
            basic.clearScreen()
            nextTask = Task.Scan // restart with a scan
            basic.showArrow(ArrowNames.West)
        break
    }

})

input.onButtonPressed(Button.B, function () {
    if (nextTask != Task.Measure) {
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
    nextTask = Task.Scan // new mode, so start with a scan
    basic.showArrow(ArrowNames.West)
})
 // normal live operation...
let config = Config.Buggy
let nextTask = Task.Scan
basic.showArrow(ArrowNames.West)
let spinRPM = 0
let turn30 = 0
