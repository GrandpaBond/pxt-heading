# pxt-heading

## Introduction
Microbit makes a great controller for a robot buggy. But every vehicle ideally needs to know which way it's pointing. 

If you know exactly how far each wheel has rotated, then from their diameter you can compute the difference 
in path-lengths. Knowing the axle-length then lets you compute changes of direction. However, only a few 
buggies have rotation sensors.

A rather simpler method is just to use a compass! The problem is that the built-in microbit compass expects to operate
when it is horizontal, yet most robot buggies mount it vertically. It also asks to be calibrated by tilting in every 
possible direction.

This pxt-heading extension is designed to meet the need for an orientation-independent compass, based just on the 
magnetometer readings. It, too, needs to first calibrate the magnetometer (to discover how it is affected by its 
magnetic environment), but that is achieved simply by spinning the buggy on the spot for a couple of rotations, 
and then teaching it where North is.

## Earth's Magnetic Field
A simplified view is that the earth's magnetic field points towards the North magnetic pole (situated incidentally 
in northern Canada!) In the southern hemisphere it points up out of the ground, and for the northern hemisphere 
it points down into the ground; near the equator it points roughly horizontally.

So, when viewed from the perspective of our spinning buggy, the magnetic field-vector sweeps out a cone:
sharply pointed in polar regions, and opened-out flat when near the equator.

## Magnetometer
Although the magnetometer provides three readings (X,Y & Z), we only need to use two to get our heading angle. 
The challenge is to choose the best two for the job! 

Depending on the specific mounting orientation of the microbit in the buggy, in a few special places and cases, 
the field-vector cone is aligned with one of the three magnetometer axes, so as the magnetic field vector 
apparently sweeps around its conical path, it traces out a neat circle on the plane of the other two axes. 
Their readings can then be easily interpreted to give the heading angle. 

However, in most places and cases this cone is tilted at an angle, and traces out an ellipse on each of the three 
orthogonal planes defined by a pair of axes (XY,YZ & ZX).  The three ellipses will in general be offset from the
origin and show different eccentricities. We will get the best heading discrimination by choosing the plane 
with the least eccentric (i.e most nearly circular) ellipse. 

Having selected the best two axes, we'll need to transform readings from around the ellipse back onto a true circle
that is centred on the origin, giving us a relative angle that can (eventually) be offset by a fixed bias to return 
the true heading with repect to North.


has a horizontal component at every location (apart from at the magnetic poles). 
As the buggy spins, this North-South component sinusoidally 
Whatever the microbit's mounting orientation, changes the flux readings in at least two of the magnetometer channels. 

These sinusoidal curves will be 90 degrees out-of-phase, but need to be corrected in various ways.

## Calibration
The first task is to determine which two of the three axes to use. Then we'll need to compensate for field-distortions 
due to the buggy itself (i.e. fixed metalwork and motor magnets close to the microbit mounting-point).

Finally we'll need to balance up the detected field-strengths (depending where you are located on the globe) so 
that we can apply simple trigonometry to compute the angular bearing with respect to North.    



## heading.scan()
It is obviously not feasible for this extension to know how to turn your buggy in any particular direction, so you 
must first set it spinning on-the-spot (by running its motors in opposite directions) and then call heading.scan(ms) 
function to perform a magnetic scan of 3-D readings. The scanning time (ms) should be sufficient to complete one or 
two full rotations.

## heading.setNorth()
Because the extension can't know just how the microbit is mounted your particular buggy, you need to point it towards 
North and then call this function to register that direction as North (or zero degrees). It can then
analyse the data collected during the scan to choose the appropriate axes and how to calibrate them. 
As an interesting by-product, it will return the spin-rate of the buggy during the previous scan (in RPM), letting 
you compare the effects of setting different motor speeds.

## heading.degrees()
Having performed the scan, and prepared the measuring set-up, this is the function that returns the current compass 
heading in degrees (0 to 360)

## rpm2Speed(diameter, axle)
A utility function to help with motor calibration. This function converts the spin-rate achieved with wheels turning 
in opposite directions, into the equivalent linear speed when both wheels are going forwards. Calculations are based 
on the wheel-diameter and axle-length. 

It should be noted that motor calibration using this technique can only ever be approximate. 
For a given power setting, inertial and frictional effects may mean that the actual 
wheel-rotation speeds achieved will differ between moving forward and spinning on the spot.
Also, for low power settings, some buggies may give an initial "kick" to get the motor going!









> Open this page at [https://grandpabond.github.io/pxt-heading/](https://grandpabond.github.io/pxt-heading/)

## Use as Extension

This repository can be added as an **extension** in MakeCode.

* open [https://makecode.microbit.org/](https://makecode.microbit.org/)
* click on **New Project**
* click on **Extensions** under the gearwheel menu
* search for **https://github.com/grandpabond/pxt-heading** and import

## Edit this project

To edit this repository in MakeCode.

* open [https://makecode.microbit.org/](https://makecode.microbit.org/)
* click on **Import** then click on **Import URL**
* paste **https://github.com/grandpabond/pxt-heading** and click import

#### Metadata (used for search, rendering)

* for PXT/microbit
<script src="https://makecode.com/gh-pages-embed.js"></script><script>makeCodeRender("{{ site.makecode.home_url }}", "{{ site.github.owner_name }}/{{ site.github.repository_name }}");</script>
