# AIRIS SIMS

## Intro
This repo is just a place to keep all the scripts used for work on AIRIS missions. You can find more information on our [website](https://washusatellite.com/), or on the page for NASA's [ADAPT mission](https://adapt.physics.wustl.edu/). Most updated script is [cmosgui2.py](#cmosgui2py) 

## Scripts
### cmossim.py

Where it all started. Now obselete, feel free to use if you aren't a fan of GUIs. Includes methods for: 
* Function to calculate limiting magnitude. Trust at your own risk.
* Random star distribution
* Gaussian PSF
* Photon Shot noise, Read Noise, Thermal Noise/Dark Current
* Clipping/Saturation Capacity
* Image Visualization (pyplot)

### cmosgui.py

Basic simulation GUI, doesn't include advanced methods, will probably deprecate soon. Includes methods for:
* Random star distribution
* Gaussian PSF
* Photon Shot noise, Read Noise, Thermal Noise/Dark Current
* Clipping/Saturation Capacity
* Image saving (png, fits)
* Image visualization (pyplot)

### cmosgui2.py

Most updated script. Will continue to update as more accurate simulation methods are used. Includes methods for:
* Everything in [cmosgui.py](#cmosguipy)
* 3x3 Pixel Binning
* Moving exposures (ie. taking long exposures while camera is moving)
* Cosmic Ray Sims (not verified for physical accuracy, only simulates track-like (straight line), not spot-like (dot) or worm-like (polyline))
* Uniform sky background (not verified for physical accuracy, see [skybackcalc.py](#skybackcalcpy) for more info)
* SNR (signal to noise ratio) calculation (uses different definitions, will update soon)

### mag_lim_test.py

Short script to try calculating magnitudes using functions from [cmossim.py](#cmossimpy). Will probably deprecate soon.

### skybackcalc.py

Attempt to calculate physically accurate signal from sky using Rayleigh scattering calculations. In its current state with current results, either extremely inaccurate or we're screwed, TBD. Will continue to update.

### snrcalc.py

Script to test SNR functions to implement in [cmosgui2.py](#cmosguipy)

## Future Goals

There's a lot of work that could potentially be done!

* Implement physically accurate sky background calculations
  * I've tried to calculate sky background signal due to atmospheric Rayleigh scattering (see skybackcalc.py) but I seem to be getting an extremely high value, one that washes out all stars. Not sure if this is unphysical or where my error is.
* Implement accurate SNR calculations
  * I've found a lot of different sources for how to accurately calculate SNR, so we need to figure out what's accurate, and what the SNR of our system is for a given magnitude. This will allow us to determine the limiting magnitude of our optical system.
* Implement calibration pipeline
  * Given these sample images, can we recover signal from the noise? We should design a calibration pipeline that takes sample images, simulated calibration frames (Flats, Darks, Bias) and the PSF map to generate cleaner images.
* Simulate accurate night sky
  * Currently, we just generate a bunch of randomly distributed point sources with randomly generated magnitudes. At some point, it may be interesting to implement some kind of global star map where we can actually model the real night sky, including some nebulae. 
  * Another addendum to this would be also simulating the exact limits of our field of view. As in, we should eventually be able to take an input direction to point in and show what the sky would look like. We're also blocked by certain objects (like the balloon we're attached to) so a lot of our FOV is occulted.

## Contact

Feel free to reach out if you have any comments/questions/recommendations or if you're just interested in the projects. 
a.wadivkar@wustl.edu
washusatellite@gmail.com
