% AKHeitman .. making things easier for future graphs
% Wraps the core computation below 


clear; clear all; close all; clc

modulation  = 'model_parameters';
mod_version = 'WN_ONFFvsLin_8pix_ID_8pix'; purpose = 'Examine effect of ONOff Hard rectification channels in WN, no cone model';
comp_params{1}.name = 'single linear channel';
comp_params{2}.name = 'hard rect ON and Off channel';
comp_params{1}.filtertype = 'fixedSP';
comp_params{2}.filtertype = 'OnOff_hardrect_fixedSP_STA';
comp_params{1}.cmodel    ='8pix_Identity_8pix';
comp_params{2}.cmodel    ='8pix_Identity_8pix';
fixed.fit_type = 'BW'; fixed.test_type = 'BW';
fixed.fitversion = 'Fit';fixed.map_type = 'mapPRJ';


cellselectiontype = 'shortlist';
exptests = [1 2 3 4];
GLMdir ='/netapp/snle/lab/Experiments/Array/Analysis/akheitman/NSEM_Projects/GLM' ;
