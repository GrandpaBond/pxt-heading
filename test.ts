// tests go here; this will not be compiled when this package is used as an extension.
input.onButtonPressed(Button.A, function () {
    switch (task) {
        case "SCAN":
            basic.showString("S")
            basic.pause(1000)
            if (!live) {
                heading.testMode(true)  // use sample data while debugging...
                heading.scan(2000)   
            } else {
                heading.testMode(false)  //
                Kitronik_Move_Motor.spin(Kitronik_Move_Motor.SpinDirections.Right, 30)
                heading.scan(4000)
                Kitronik_Move_Motor.stop()
            }
            basic.showIcon(IconNames.Yes)
            task = "ANALYSE"
        break

        case "ANALYSE":
            basic.showString("A")
            basic.pause(1000)
            basic.clearScreen()
            spinRPM = heading.analyseScan()
            basic.showNumber(spinRPM)
            turn30 = 60000 / (12 * spinRPM) // time to turn 30 degree
            datalogger.setColumnTitles("uRaw", "vRaw", "u", "v", "val")
            task = "MEASURE"
            basic.pause(3000)
            basic.showIcon(IconNames.Yes)
            basic.pause(500)
            basic.showString("?")
        break

        case "MEASURE":
            basic.showIcon(IconNames.No)
            basic.pause(1000)
            basic.clearScreen()
            task = "SCAN" // restart with a scan
        break
    }

})

input.onButtonPressed(Button.B, function () {
    if (live){
        Kitronik_Move_Motor.spin(Kitronik_Move_Motor.SpinDirections.Right, 30)
        basic.pause(turn30)
        Kitronik_Move_Motor.stop()
    }
    let compass = heading.degrees()
    basic.clearScreen()   
    basic.pause(500)
    basic.showNumber(Math.floor(compass))
    basic.pause(1000)
    basic.clearScreen()   
    basic.pause(500)
    basic.showString("?")
})


input.onLogoEvent(TouchButtonEvent.Pressed, function () {
    if (live) {
        live = false
        basic.showIcon(IconNames.Skull)
    } else {
        live = true
        basic.showArrow(ArrowNames.North)
    }
    basic.pause(3000)
    basic.clearScreen()
    task = "SCAN"
})

let live = true
let task = "SCAN"
basic.showArrow(ArrowNames.North)
basic.pause(3000)
basic.clearScreen()
let spinRPM = 0
let turn30 = 0
