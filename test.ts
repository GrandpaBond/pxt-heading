// tests go here; this will not be compiled when this package is used as an extension.
input.onButtonPressed(Button.A, function () {
    basic.showIcon(IconNames.Heart)
    basic.pause(1000)
    if (testing) {
        heading.simulateScan()  // use sample data while debugging...   
    } else {
        Kitronik_Move_Motor.spin(Kitronik_Move_Motor.SpinDirections.Right, 30)
        heading.scan(4000)
        Kitronik_Move_Motor.stop()
        // heading.dumpData()
        basic.clearScreen()
    }
})
input.onButtonPressed(Button.B, function () {
    basic.showString("?")
    spinRPM = heading.prepare(true)
    //heading.dumpLimits()
    basic.clearScreen()
    basic.showNumber(spinRPM)
    datalogger.setColumnTitles("uRaw", "vRaw", "u", "v", "val")
    basic.pause(1000)
})


input.onLogoEvent(TouchButtonEvent.Pressed, function() {
    let turn30 = 60000 / (12 * spinRPM) // time to turn 30 degrees
    for (let i = 0; i < 24 ; i ++) {  // twice round the clock
        Kitronik_Move_Motor.spin(Kitronik_Move_Motor.SpinDirections.Right, 30)
        basic.pause(turn30)
        Kitronik_Move_Motor.stop()
        let compass = heading.degrees()
        basic.showNumber(Math.floor(compass))
        basic.pause(200)
    }
})

let testing = true
let test = 0
let spinRPM = 0
basic.pause(1000)
basic.clearScreen()
