---------------------------------
SPIKE SORTING README 
(tamachado 2009-06-23):
---------------------------------

This project folder contains the following code/resources about how to use MATLAB in order to interface with each step of the spike sorting algorithm. In this document are references to information about the current spike sorting framework, ideas for improvements, and descriptions of the code stored in this folder.

If any necessary MATLAB-Vision interfaces are written and added to Vision, they should be new classes within the matlab package (edu.ucsc.neurobiology.vision.matlab).

---------------------------------
1. Spike Sorting/Vision Interface
---------------------------------

These are the principal steps of the current spike sorting technique along with the interfaces in MATLAB:

spike finding: raw data reading
PCA/projections generation: projections file reading
clustering: /imm-util/, model file reading/writing
neuron cleaning: /neuron-cleaning/, neuron file reading/writing

Each of these interfaces with Java are heavily documented in ../documentation/common_java_calls.m

That file contains a series of examples for how to read and write to nearly every intermediate file generated during spike sorting. 

YOU SHOULD LOOK AT THAT FILE BECAUSE IT HAS EXAMPLES SHOWING HOW TO ACCOMPLISH MOST TASKS THAT YOU MIGHT WANT TO DO.

---------------------------------
2. Code Stored in This Directory
---------------------------------

imm-cluster/ -- a fully functioning framework (along with a script, imm-auto-cluster) that runs an advanced (a variant of Rasmussen's Infinite Gaussian Mixture Model written by Frank Wood) clustering algorithm on a given dataset. This framework is pretty well-written and might be useful as a starting point for other spike sorting projects. If someone wants to write a new clustering algorithm in MATLAB, this framework already does basically everything that you could want to do. The script to run it is already on the server in the scripts folder, so just type 'imm-auto-cluster' or examine the source in the subfolder of this directory.

duplicate-removal/ -- some unfinished work was done to try and improve duplicate removal. This code, along with some figures about the state of the project, are present in this directory. While the algorithms were not completely tuned up, the code should work fine and is definitely a good starting point for further exploration into duplicate finding.

references/ -- contains a (slightly outdated) document about the state of the spike sorting algorithm that documents various attempts at improving things that have been tried. It also contains a copy of Maneesh Sahani's thesis, which has some good, straightforward ideas for improving spike sorting (such as ellipsoidal thresholding and noise whitening). A copy of Frank Wood's paper about his spike sorting/clustering algorithm is also included.

---------------------------------
3. Next Steps/Project Ideas
---------------------------------
Assuming the current spike sorting framework is improved (instead of replaced) over the next few years, here is a list of improvements that could be implemented reasonably quickly--and that should be tried at some point:

1. Improve duplicate removal: consider merging duplicate cells instead of throwing them away. While duplicate removal is a reasonably difficult problem, it is definitely possible to find many more spikes per cell if just the unique spikes present in each member of clear duplicate pairs are kept (instead of just keeping one of the two duplicate copies). See work in the duplicate-removal subfolder.

2. Iterative PCA. By performing PCA multiple times on subsets of the data, it is possible to extract many extra cells. The downside of this serial technique is that it is order dependent. However, by automating this process, and tuning the parameters that determine how the algorithm chooses  subsets, it's likely that work implementing such an algorithm would be useful.

3. Improving thresholding. The way vision currently applies thresholds to the raw data during spike finding is far from optimal. There are certainly better ways to do thresholding, such as thresholding across electrodes (ie the threshold on electrode k is the sum of the squared amplitudes on electrode k and its neighbor electrodes). This is described in Maneesh Sahani's thesis (which is included in the references folder). 

4. Resolving Superposition Spikes. There currently exists a framework developed by Jonathan Pillow and Jon Shlens to find and extract spikes, even when multiple cells were firing simultaneously. However, this is a pretty slow process, and it's likely that some efficient (and likely simplified) version of this analysis could be implemented for day to day use. This would probably be a pretty hard and time consuming project, but it's one that should be done at some point.


