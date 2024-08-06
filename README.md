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
based solely on the magnetometer readings. It does not (yet) include tilt-compensation, so it expects the buggy to be operating on a flat surface.

The magnetometer does still need to be calibrated (to discover how it 
sees its magnetic environment), but we achieve that simply by spinning the buggy on the spot for a few rotations, 
and then teaching it where North is.



After the block definitions below, you can read more about The Earth's Magnetic Field; Magnetometer Calibration; and Field Distortions.

## Scanning Around

```sig
heading.scanClockwise(ms): number
```
Before navigating by compass headings, your buggy must first perform a magnetic scan of its environment.
It is obviously not feasible for this extension to know how to turn every buggy in any particular direction,
so you must first set it spinning clockwise on-the-spot (by running its motors in opposite directions),
and then call the ``||heading:scanClockwise()||`` function to take 3-D magnetometer readings all around.
When enough data has been collected it will be analysed, and if OK, calibration parameters will be set up
for taking future readings.

> ``||heading:ms||`` is the scanning time in ms (sufficient to complete about three full rotations).

Returns zero if successful, or a negative error code:

- -1 : NOT ENOUGH SCAN DATA

- -2 : FIELD STRENGTH TOO WEAK

- -3 : NOT ENOUGH SCAN ROTATION


### ~reminder
    If your buggy turns in small circles rather than spinning on-the-spot, you'll need to balance the power of the two 
    wheel motors by adjusting its motor bias.
### ~


## Where's North?

```sig
heading.setNorth()
```

Because this extension can't know exactly how the microbit is mounted on your particular buggy, you'll need to physically point it towards 
"North" and then call this function (after first performing a ``||heading:scanClockwise()||`` ).
It will then take a fix on the current heading and register that direction as zero degrees.

### ~reminder
    The actual direction of the buggy when this function is called is arbitrary: it could be Magnetic North; or True North (compensating for local declination); or any convenient direction from which to measure subsequent heading angles.
### ~

## Where am I pointing?

```sig
heading.degrees(): number
```

Having first calibrated the measuring set-up using ``||heading:scanClockwise()||``,
and then specified the zero-degrees direction using ``||heading:setNorth()||``, 
this is the function that actually returns the current compass heading in degrees (0 to 360), clockwise from "North".

### ~reminder
    NOTE: Readings are only accurate to a few degrees, and are very sensitive to variations in the tilt of the buggy!
### ~

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
heading.equivalentSpeed(axleLength: number): number
```

A utility function to help (somewhat) with motor calibration, this function converts the spin-rate achieved with wheels turning 
in opposite directions into the equivalent linear speed (in mm/sec) if both wheels were going forwards. 
To calculate this we just need to know how far apart the wheels are.

> ``||heading:axleLength||`` is the wheel separation (in mm). 

### ~reminder
It should be noted that motor-calibration using this technique can only ever be approximate. 

For a given power setting, inertial effects and tyre-friction may mean that the wheel-rotation speeds
actually achieved will differ between moving forwards and spinning on the spot.

Also, on slower power settings, some buggies may complicate matters by automatically giving an initial high-power "kick", 
just to get the motors going!
### ~

# Background Information
This section delves a bit deeper into the physics and maths underpinning this extension.

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
The microBit magnetometer reads the field strength in three independent directions (X, Y & Z), giving a 3-D vector. 
You could think of the three readings as the dimensions of a 3D box with the local magnetic field-vector as its diagonal. 
Some part of this measurement will be due to Earth's magnetic field, though there may be contributions from other magnetic sources 
(see "Field Distortions "below)

In our calibration approach, we imagine the situation from the perspective of a buggy spinning on the spot.
As the buggy turns clockwise, the end of the magnetic field-vector appears to sweep out an anti-clockwise path we call the Spin-Circle.

Although the magnetometer provides three readings (X,Y & Z), we will only ever be needing to use two of them to fix our 
heading angle; the challenge is to choose the best two for the job! 

### Square Mounting
If the microbit is mounted squarely in the buggy, the axis of the Spin-Circle will be neatly aligned with whichever
magnetometer axis is pointing "up". The simple sine-wave readings taken from the remaining two axes can then be easily 
interpreted to give us the heading angle, remembering that the magnetometer chip is mounted on the **back** of the microbit, 
so its axes are effectively reversed!

### Slanted Mounting
However, in the fully general situation where the microbit might be mounted on a slant (e.g. as a "dashboard"), the Spin-Circle 
will appear fore-shortened into an Ellipse when looking down (or up) any particular axis onto the plane of the other two axes (YZ, ZX or XY).

This Ellipse will have an arbitrary centre, a degree of eccentricity, and maybe a tilt. Fore-shortening means that
equally-spaced heading angles will appear bunched around the pointed ends of the Ellipse. 

To recover the buggy's correct heading, we'll need to undo this, mapping points on the Ellipse back onto the Spin-Circle.

Using the calibration data, we compute the different Ellipse eccentricities in each of the three possible views (YZ, ZX or XY).
We will always get the most accurate heading discrimination from the "roundest", most "square-on" view (i.e. the one with the 
least eccentric Ellipse).

### Transforming Readings
Having selected the best view to use (and thus the two most useful axes), we'll need to transform each new 2-D reading from its
fore-shortened position on the Ellipse to its equivalent position on the Spin-Circle. This can require up to four transformation steps:

1) First we will need to shift each coordinate by a fixed offset (effectively moving the Ellipse's centre to the Origin), 

2) If we are unlucky and the Ellipse appears tilted, we may then have to rotate our reading so that the Ellipse lies squarely 
over the axes;

3) Then we must scale one value in proportion to the eccentricity (so stretching the Ellipse back into a perfect circle);

4) Finally, for a tilted Ellipse, we must undo the axis-alignment rotation in step 2. 

Taking the the angle of this final "stretched" position on the circle and subtracting the "North" angle then gives us the true heading.

### Finding Ellipse Properties

This correction can in theory be applied in any of the three views (unless exactly side on to
the Spin-Circle), but the best accuracy is obtained from the most circular of the three views.
Readings on a near-circular Ellipse are barely fore-shortened at all, so we can skip correction!

So for each view we must derive the two important Ellipse properties: {tilt} and {eccentricity}. 
This first requires detection of its major and minor axes. The maths for fitting an ellipse to noisy
2D data is both complex and fairly inaccurate. Luckily we have a sequence of 3D readings and the 
third orthogonal dimension (the Normal) lets us use a simpler analysis method to find the axes.

So for example: in the XY plane, the ellipse radius {x^2 + y^2} will be at a maximum as it is passing its major-axis. 
This happens twice per rotation, when the field is crossing the XY plane. 

We know that the three magnetometer readings are related by the formula:  {x^2 + y^2 + z^2 = B}, where B is the 
constant (over our time-scales) magnetic field, so near each end of the major-axis the z-value will be at a minimum:
at one end it dips below the plane; at the other end it rises above.

Conversely the XY radius is at a minimum (its minor-axis) where the field points farthest from the plane.
Here the z-value is at its biggest, and its sign will differ between the two ends of the XY minor-axis.

The same holds true for the other two planes: the Normal helps us sorts out which end is which for each axis. 
For each of these three mutually-orthogonal views, we re-label its coordinates as {u,v}, with {w} being the
third (orthogonal) coordinate.


## Field Distortions
So much for the theory! In the real world there are several factors which conspire to distort the actual local field, as 
measured by the magnetometer. These can include:

* Departure from horizontal: any tilting of the buggy due to uneven ground will seriously distort the readings.

* Environmental magnetic anomalies, due to magnets or electrical machinery near where the buggy is being used.

* Electro-magnetic fields due to flowing currents (including those powering the buggy's motors or servos).

* The magnetometer chip is actually temperature-sensitive, so readings may drift as it gets warmer with repeated use.

* Varying fields from the arbitrary angles of any rotating magnetic elements in the motors at the moment of measurement.

* Fixed metalwork and static motor or loudspeaker magnets on the buggy itself that are positioned close to the microbit. 
These so-called "hard-iron" distortions will add (or subtract) a different fixed bias for each axis.

* "Jitter". Even when nothing has changed, repeated readings from the magnetometer often tend to differ.

We can only ever hope to account for the last two of these errors. Spinning the buggy lets us separate the local and environmental
influences and so calculate the fixed biases to be applied to each axis reading. Additionally, to smooth out any jitter,
we always take several consecutive readings and average them.

Counter-intuitively, reading the heading while the motors are running may actually prove **more** accurate, as their rotating fields 
could then get averaged-out somewhat.

Unfortunately, many of these factors lie outside our control, and will inevitably limit the accuracy and repeatability
of reported headings to within five degrees or so.


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