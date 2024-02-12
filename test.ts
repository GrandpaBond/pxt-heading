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
    heading.dumpLimits()
    basic.clearScreen()
    basic.showNumber(spinRPM)
    basic.pause(1000)
})


input.onButtonPressed(Button.AB, function () {
    basic.showIcon(IconNames.Heart)
    //meter.use(meter.Styles.Dial, 0, 360)
    while (test < 300) {
        let compass = heading.degrees()
        //basic.showNumber(Math.floor(compass))

        //meter.show(compass,250)
    }
})

let testing = true
let test = -1
basic.pause(1000)
basic.clearScreen()
