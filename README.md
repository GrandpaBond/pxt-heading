```package
heading=github:grandpabond/pxt-heading
```

# Heading Extension

## Introduction
MicroBit makes a great controller for a robot buggy. But every vehicle ideally needs to know which way it's pointing. 

If you know exactly how far each wheel has rotated, then from their diameter you can compute the difference 
in path-lengths. Knowing the axle-length then lets you compute changes of direction. However, only a few 
buggies have rotation sensors...

A rather simpler method is just to use the compass! The problem is that the built-in microbit compass is fairly imprecise 
and seems to operate best when it is horizontal, yet most robot buggies mount it vertically. 
It also demands calibration by tilting it in every possible direction --potentially a rather awkward operation with a buggy!

This ``||heading:heading||`` extension is designed to meet the need for a location-independent and orientation-independent compass,
based solely on the magnetometer readings. The magnetometer does still need to be calibrated (to discover how it 
sees its magnetic environment), but we achieve that simply by spinning the buggy on the spot for a few rotations, 
and then teaching it where North is.

## Earth's Magnetic Field
A simplified view is that the Earth's magnetic field loops round from the magnetic South pole towards the magnetic North 
pole. Because the field is generated by circulation currents in the Earth's liquid iron core, the situation is actually much more complicated, 
and these poles are gradually wandering.

A decade ago, the Magnetic north pole was situated in northern Canada, some 500km from the true pole; it is currently heading for 
Siberia at around 50 km per year! So in reality a compass seldom points North, but it can however be relied on to change only slowly
over time.

In the southern hemisphere the field-vector points up out of the ground, and for the northern hemisphere it points down 
into the ground; near the Equator it is roughly horizontal. 


## Magnetometer Calibration
The microBit magnetometer reads the field strength in three independent directions (X, Y & Z), giving the dimensions of a 3D box 
with the local magnetic field-vector as its diagonal. Some part of this measurement will be due to Earth's magnetic field.

In our calibration approach, we imagine the situation from the perspective of a buggy spinning on the spot.
As the buggy turns clockwise, the end of the magnetic field-vector sweeps out an anti-clockwise path we call the Spin-Circle.

Although the magnetometer provides three readings (X,Y & Z), we will only ever be needing to use two of them to get our 
heading angle; the challenge is to choose the best two for the job! 

### Square Mounting
If the microbit is mounted squarely in the buggy, the axis of the Spin-Circle will be neatly aligned with whichever
magnetometer axis is pointing "up". The simple sine-wave readings taken from the remaining two axes can then be easily interpreted to give 
us the heading angle. 

### Slanted Mounting
However, in the fully general situation where the microbit might be mounted on a slant, the Spin-Circle appears 
fore-shortened into an Ellipse when looking down any particular axis onto the plane of the other two axes (YZ, ZX or XY).

This Ellipse will have an arbitrary centre, a degree of eccentricity, and maybe a tilt. Fore-shortening means that
equally-spaced heading angles will appear bunched around the pointed ends of the Ellipse. 

To recover the buggy's correct heading, we'll need to undo this, mapping points on the Ellipse back onto the Spin-Circle.

Using the calibration data, we compute the different Ellipse eccentricities in each of the three possible views (YZ, ZX or XY).
We will always get the most accurate heading discrimination from the "roundest", most "square-on" view (i.e. the one with the 
least eccentric Ellipse).


### Transforming Readings
Having selected the best view to use (and thus the two most useful axes), we'll need to transform each new 2-D reading from its
fore-shortened position on the Ellipse to its equivalent position on the Spin-Circle. This can require up to four transformation steps:


1) First we will need to shift each coordinate by a fixed offset (effectively moving the Ellipse's centre to the origin), 

2) If we are unlucky and the Ellipse appears tilted, we may then have to rotate our reading so that the Ellipse lies squarely 
over the axes;

3) Then we must scale one value in proportion to the eccentricity (so stretching the Ellipse back into a perfect circle);

4) Finally, for a tilted Ellipse, we must undo the axis-alignment rotation in step 2 

Taking the the angle of this final "stretched" position on the circle and subtracting the "North" angle then gives us the true heading.

## Field Distortions
So much for the theory! In the real world there are several factors which conspire to distort the actual local field, as 
measured by the magnetometer. These can include:

* Environmental magnetic anomalies, due to magnets or electrical machinery near where the buggy is being used.

* Electro-magnetic fields due to flowing currents (including those powering the buggy's motors or servos).

* Varying fields from the arbitrary angles of any rotating magnetic elements in the motors at the moment of measurement.

* Fixed metalwork and static motor or loudspeaker magnets on the buggy itself, if they lie close to the microbit mounting-point. 
These so-called "hard-iron" distortions will add (or subtract) a different fixed bias for each axis.

* "Jitter". Even when nothing has changed, repeated readings from the magnetometer often tend to differ.

We can only ever hope to account for the last two of these errors. Spinning the buggy lets us separate the local and environmental
influences and so calculate the fixed biases to be applied to each axis reading. Additionally, to smooth out any jitter,
we always use an average of seven consecutive readings.

Counter-intuitively, reading the heading while the motors are running may actually prove **more** accurate, as their rotating fields 
will then get averaged-out somewhat.

Unfortunately, the first two factors lie outside our control, and will inevitably limit the accuracy and repeatability
of reported headings.


## Scanning Around

```sig
heading.scanClockwise(ms)
```

It is obviously not feasible for this extension to know how to turn every buggy in any particular direction, so you 
must first set it spinning clockwise on-the-spot (by running its motors in opposite directions)
and then call the ``||heading:scanClockwise()||`` function to perform a magnetic scan of 3-D readings. 

> ``||heading:ms||`` the scanning time in ms (sufficient to complete at least two full rotations).

## Where's North?

```sig
heading.setNorth(): number
```

Because this extension can't know just how the microbit is mounted on your particular buggy, you'll need to physically point it towards 
"North" and then call this function (after first performing a ``||heading:scanClockwise()||`` ).
This function will first analyse the data collected during the previous scan (to choose the best axes and how to calibrate 
readings taken from them).
It will then take a fix on the current heading and register that direction as zero degrees.

## Where am I pointing?

```sig
heading.degrees(): number
```

Having performed the scan using ``||heading:scanClockwise()||``, and calibrated the measuring set-up using ``||heading:setNorth()||`` , 
this is the function that actually returns the current compass heading in degrees (0 to 360) clockwise from North.

## Rotation Time

```sig
heading.spinTime(): number
```

As an interesting by-product of analysing a scan, we calculate the rotation-period of the buggy. This function returns
that time (in ms), so letting you compare the effects of scanning with different motor speed settings.

## Rotation Rate

```sig
heading.spinRate(): number
```

This function returns the latest scan rotation rate, expressed in revs-per-minute (RPM)

## Power v. Speed

```sig
equivalentSpeed(axleLength: number): number
```

A utility function to help (a little) with motor calibration, this function converts the spin-rate achieved with wheels turning 
in opposite directions into the equivalent linear speed (in mm/sec) if both wheels were going forwards. 

> ``||heading:axleLength||`` calculations require knowledge of the wheel separation (in mm). 
