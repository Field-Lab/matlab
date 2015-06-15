Current version (06/15/2015) is organized as follows:
Directory
/Users/gomena/Research/GIT/ChichilniskyGIT/matlab/code/projects/electrical_stim/spike_sorting_GM/
in the main folder there are the following files
-ExampleSeveralPatterns.m
-ExampleSeveralElectrodes.m
-Documentation.rtf
The example files are self contained and show how to run the algorithm, with comments in all relevant steps. The documentation.rtf Is important to understand the details of examples and the functions.

Also, there are folders
Utils: Auxiliary functions
Core: The main functions of the algorithms
Display: Code for displaying results
ResultsAnalysis: Code for analizing results (in principle, from the 06/07/2015 week with Lauren and Gonzalo in Stanford)
NumericalExperimentsGonzalo

Please don't touch the core/utils functions before discussing them. There are still a lot of uncommented functions, I plan on creating help for them soon
Gonzalo (06/15/2015)


Gonzalo (06/15/2015)
New function SpikeSortingCompact was created, see the help for usage. The aim is to don't have to run examples, but do Spike Sorting with just having to fill the more esentiall variables, the path and neurons, and in such a way that all values are specified as defaults (as stated in the examples and fillDefaultValues function) unless otherwise stated. Also, no need to translate templates to have them aligned to onset of spike at time 11. (using the also new AlignTemplates function)
This function should be easily modified in order to give more flexibility to define default values. 
