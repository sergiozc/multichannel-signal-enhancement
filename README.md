# multichannel-signal-ehancement

An implementation of a MVDR beamformer for a linear array of acoustic sensors can be found at this repository. The corresponding code is found at "./macros/main.m".

The target signal is a real signal for which the sound source is 1 meter from the array.
The spherical wave and plane wave approximation have been compared. Moreover, the MVDR beamformer has been compared with the Delay & Sum beamformer.

Finally, an additional version has been implemented in which different sub-array configurations are selected with the purpose of optimizing the output as a function of frequency. This version can be found at "./macros/version_subarrays.m".

All the results (.wav format) are within "RESULTS" directory.
