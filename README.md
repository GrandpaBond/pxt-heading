# pxt-heading

## Introduction
Microbit makes a great controller for a robot buggy. But every vehicle ideally needs to know which way it's pointing. 

If you know exactly how far each wheel has rotated, then from their diameter you can compute the difference 
in path-lengths. Knowing the axle-length then lets you compute changes of direction. However, only a few 
buggies have rotation sensors.

A rather simpler method is just to use a compass! The problem is that the built-in microbit compass expects to operate
when it is horizontal, yet most robot buggies mount it vertically. It also asks to be calibrated by tilting it in every 
possible direction, potentially a rather awkward operation with a buggy...

This pxt-heading extension is designed to meet the need for a location-independent and orientation-independent compass,
based solely on the magnetometer readings. 
The magnetometer does still need to be calibrated (to discover how it 
sees its magnetic environment), but we achieve that by simply spinning the buggy on the spot for a couple of rotations, 
and then teaching it where North is.


## Earth's Magnetic Field
A simplified view is that the earth's magnetic field points from the South magnetic pole towards the North magnetic 
pole (situated incidentally in northern Canada!) In the southern hemisphere it points up out of the ground, 
and for the northern hemisphere it points down into the ground; near the equator it points roughly horizontally.

The magnetometer reads the field strength separately in the X, Y & Z directions (the dimensions of a 3D box with the magnetic field vector as its diagonal). As the buggy spins,these orthogonal strengths oscillate either side of zero.


One way of visualising the situation is to imagine a sphere centred on the magnetometer. A magnetic field-vector 
through the middle will intersect it at a point on its surface. When viewed from the perspective of our spinning buggy, successive points will traces out a circle on the sphere whose diameter ranges from tiny (in polar regions where the field is
more vertical) to the full diameter of the sphere (when we are near the earth's equator and the field is horizontal).
We'll call this circle the Spin Circle.

## Magnetometer Calibration

Although the magnetometer provides three readings (X,Y & Z), we will only need to use two of them to get our 
heading angle; the challenge is to choose the best two for the job! 

Depending on the specific mounting orientation of the microbit in the buggy, we may find that the axis of the 
Spin Circle will be neatly aligned with one of the three magnetometer axes. The readings of the remaining two 
axes can then be easily interpreted to give us the heading angle. 

However, in the fully general situation where the microbit might be mounted on a slant, the Spin Circle appears 
tilted into an ellipse when looking down one or more axes onto the plane of the other two axes (XY,YZ & ZX)
Each ellipse will in general show a different eccentricity and we will obtain the best heading discrimination by 
choosing the projection plane showing the least eccentric (i.e most nearly circular) ellipse. 

(Some systems just use an accelerometer to tell which way is "up", but this can only be trusted if the buggy is 
guaranteed to be completely at rest when you measure.)

and also (typically) be offset from the origin

Having selected these best two axes, we'll need to transform each new reading from around the ellipse back onto 
a true circle. This can require up to three steps:

* first we will generally need to shift it, effectively moving the ellipse's centre to the origin; 

* if we are unlucky and the ellipse is tilted, we may then have to rotate it to lie squarely over the axes;

* finally we must scale one value in propotion to the eccentricity (stretching the ellipse into a perfect circle).

Applying simple trigonometry to the final position on the circle gives us a relative angle that can be offset 
by a known bias to return the true heading with repect to North.

## Field Distortions
So much for the theory! In the real world there are several factors which conspire to distort the actual field as 
measured by the magnetometer:

* Environmental magnetic anomalies, due to magnets or electrical machinery near where the buggy is being used.

* Electro-magnetic fields due to flowing currents (especially to the buggy's motors).

* Fixed metalwork and static motor magnets on the buggy itself close to the microbit mounting-point. These will
    add (or subtract) a different fixed bias for each axis.

* "Jitter". Even when nothing has changed, repeated readings from the magnetometer often tend to differ.

We can only really account for the last two of these: Spinning the buggy lets us separate the local and environmental
influences and so calculate the fixed biases to be applied to each axis reading. Additionally, to smooth out any jitter,
we always use a rolling average of seven consecutive readings.

Unfortunately, the other factors sre outside our control, and will inevitably limit the accuracy and repeatability
of reported headings.

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
Having performed the scan using heading.scan(), and calibrated the measuring set-up using heading.setNorth(), 
this is the function that returns the current compass heading in degrees (0 to 360)

## rpm2Speed(diameter, axle)
A utility function to help with motor calibration. This function converts the spin-rate achieved with wheels turning 
in opposite directions, into the equivalent linear speed when both wheels are going forwards. Calculations are based 
on the wheel-diameter and axle-length. 

NOTE: It should be noted that motor calibration using this technique can only ever be approximate. 
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