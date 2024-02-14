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
    let spinRPM = heading.prepare(true)
    //heading.dumpLimits()
    basic.clearScreen()
    basic.showNumber(spinRPM)
    basic.pause(1000)
})


input.onButtonPressed(Button.AB, function () {
    basic.showIcon(IconNames.Heart)
    //meter.use(meter.Styles.Dial, 0, 360)
    while (test < 300) {
        let compass = heading.degrees()
        basic.showNumber(Math.floor(compass))
        basic.pause(1000)
        //meter.show(compass,250)
    }
})

input.onLogoEvent(TouchButtonEvent.Pressed, function() {
    let compass = heading.degrees()
    basic.showNumber(Math.floor(compass))
})

let testing = true
let turning = true
let test = 0
basic.pause(1000)
basic.clearScreen()
