MEA Data Processing Software
==============================
This is a set of scripts meant for processing of multi-electrode array (MEA)
data. The scripts interface with Vision and some will required access to the
Vision jar archive in order to work. They provide the following functionality:
   * Electric artifact estimation and removal
   * Neuron firing statistics computation
   * Neuron firing statistics plotting
   * Interface with the data acquisition software
   * Interface with Vision
   * Positioning utilities to switch between array, implant and white noise (WN)
   coordinates

Artifact estimation and removal
------------------------------
The [Artifact](./artifact) module provides functionality for estimating and 
removing from raw data electric stimulation artifacts.

Neuron statistics computation
------------------------------
The statistics computation functions may be found in the module [Statistics](./statistics).
The following statistics can be computed:
   * Peristimulus time histograms.
   * Activation data (number of spikes elicited and suppressed).
   * Electric receptive fields.
   * Estimation of the response to a spot of visible light from WN data.
   * Estimation of non-linear summation effects in the case of double spot
   stimulation.

Neuron statistics plotting
------------------------------
The plotting functions may be found in the module [Plotting](./plotting).
The following plots can be obtained:
   * Peristimulus time histograms.
   * Sigmoid activation curves.
   * Still and time-dependent electric receptive fields.
   * Neuron activation maps over the array.

Interface with the data acquisition software
------------------------------
The module [Daq Interface](./daq_interface) provides functionality for 
interfacing with the stimulation software logfiles.

Interface with Vision
------------------------------
The module [Vision Interface](./vision_interface) provides functionality for 
interfacing with Vision, by reading, writing vision files, or directly pooling
data and processing some of them with for example the removeDuplicates() 
function, which replaces some of Vision's functionality.

Other modules
------------------------------
[Resources](./resources) has  information about array and implant electrode 
layouts, which are necessary in order to map coordinates between WN, array and
implant data sets.

[Mapping](./mapping) provides the functionality for switching from one 
coordinate system to another.

[Scripts](./scripts) provides examples on how to call all the functions.