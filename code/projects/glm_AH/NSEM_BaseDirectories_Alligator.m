% AKHeitman 2014-04-02

%%%%% CAUTION CHNAGES HER WILL BE PROPAGATED THROUGHOUT THE CODE %%%%%
%%%% NEEDS TO GET UPDATED EACH TIME WE MOVE %%%%%
% This will get called throughout at will.. All changes here

% Define the important base directories that will get called
% The output saving should follow directly from there

% Try to update each time we move things around
% Try to update each time we generate new directory


function [BASEDIR] = NSEM_BaseDirectories_Alligator

%javaaddpath('/Users/akheitman/Dropbox/Lab/Applications/Vision.app/Contents/Resources/Java/Vision.jar');

BASEDIR.note                = 'GLM and NSEM directories.  Code and Data.  Independent of experiment number, parameters, fit types';
BASEDIR.NSEM_home_snleapp   = '/netapp/snle/lab/Experiments/Array/Analysis/akheitman/NSEM_Projects';

BASEDIR.general_codehome = '/Users/akheitman/Matlab_code/';
BASEDIR.analysisdir      = '/Users/akheitman/temp_Analysis';
BASEDIR.NSEM_home       = '/Users/akheitman/NSEM_Home';




BASEDIR.BlockedSpikes   = '/Users/akheitman/NSEM_Home/BlockedSpikes';
BASEDIR.Cell_Selection  = '/Users/akheitman/NSEM_Home/Cell_Selection';
BASEDIR.GLM_output_raw  = '/Users/akheitman/NSEM_Home/GLM_Output_Raw';
BASEDIR.GLM_output_analysis  = '/Users/akheitman/NSEM_Home/GLM_Output_Analysis';
BASEDIR.GLM_codehome    = '/Users/akheitman/Matlab_code/glm_AH/';





BASEDIR.GLM_troubleshootplots = sprintf('%s/troubleshootingplots', BASEDIR.GLM_codehome); 


BASEDIR.tempnote = 'temporary directory discrepancies due to Stanford move' ;
BASEDIR.temp.GLM_output_netappsnle = '/netapp/snle/lab/Experiments/Array/Analysis/akheitman/NSEM_Projects/GLM'; 
BASEDIR.temp.GLM_output_alligator = '/Users/akheitman/Matlab_code/glm_output';

% Directories that should follow over
BASEDIR.NSEM_stimuli = sprintf('%s/Stimuli', BASEDIR.NSEM_home);

end


