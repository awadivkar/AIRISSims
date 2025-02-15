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

## Contact

Feel free to reach out if you have any comments/questions/recommendations or if you're just interested in the projects. 
a.wadivkar@wustl.edu
washusatellite@gmail.com
