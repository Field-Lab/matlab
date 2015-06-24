% NBrackbill 2015-04-20

%%%%% CAUTION CHNAGES HER WILL BE PROPAGATED THROUGHOUT THE CODE %%%%%
%%%% NEEDS TO GET UPDATED EACH TIME WE MOVE %%%%%
% This will get called throughout at will.. All changes here

% Define the important base directories that will get called
% The output saving should follow directly from there

% Try to update each time we move things around
% Try to update each time we generate new directory


function [BASEDIR] = NSEM_BaseDirectories_Lovelight


BASEDIR.note                = 'GLM and NSEM directories.  Code and Data.  Independent of experiment number, parameters, fit types';
BASEDIR.NSEM_home_snleapp   = '/netapp/snle/lab/Experiments/Array/Analysis/akheitman/NSEM_Projects';

%BASEDIR.general_codehome = '/Users/Nora/Documents/MATLAB/matlab/code/';
BASEDIR.general_codehome = '/home/vision/Nora/matlab/code/'; % BERTHA
BASEDIR.analysisdir      = '/Volumes/Analysis';
BASEDIR.NSEM_home        = '/Volumes/Lab/Users/akheitman/NSEM_Home';

BASEDIR.BlockedSpikes         = '/Volumes/Lab/Users/akheitman/NSEM_Home/BlockedSpikes';
BASEDIR.Cell_Selection        = '/Volumes/Lab/Users/akheitman/NSEM_Home/Cell_Selection';
BASEDIR.GLM_output_raw        = '/Volumes/Lab/Users/Nora/NSEM_Home/GLMOutput_Raw';
BASEDIR.GLM_output_analysis   = '/Volumes/Lab/Users/Nora/NSEM_Home/GLMOutput_Analysis';


BASEDIR.GLM_develop_output_raw        = '/Volumes/Lab/Users/Nora/NSEM_Home/GLM_Develop_Output';

BASEDIR.GLM_codehome          = [BASEDIR.general_codehome 'projects/glm/'];
BASEDIR.GLM_troubleshootplots = sprintf('%s/troubleshootingplots', BASEDIR.GLM_codehome); 

% 
% BASEDIR.tempnote = 'temporary directory discrepancies due to Stanford move' ;
% BASEDIR.temp.GLM_output_netappsnle = '/netapp/snle/lab/Experiments/Array/Analysis/akheitman/NSEM_Projects/GLM'; 
% BASEDIR.temp.GLM_output_alligator = '/Users/akheitman/Matlab_code/glm_output';

% Directories that should follow over
BASEDIR.NSEM_stimuli = sprintf('%s/Stimuli', BASEDIR.NSEM_home);

end


