% Java library
javaaddpath('/Volumes/Lab/Development/vision7/Vision.app/Contents/Resources/Java/Vision.jar');

% matlab code paths
addpath(genpath('/home/vision/Nora/matlab/code'));
addpath(genpath('/home/vision/Nora/matlab/code/projects/glm'));

% set some default plot stuff
set(0, 'DefaultAxesFontSize', 18, 'DefaultAxesFontName', 'Helvetica')

datarun = load_data('2005-04-26-1/data005-pca/data005/data005');
datarun = load_params(datarun);
on_par = datarun.vision.cell_types{1,33}.cell_ids;

on_par2 = on_par(on_par > 5626);


fittedGLM = glm_fit_from_WN(on_par2, '2005-04-26-1/data005-pca/data005/data005', 'BW-20-1-0.48-11111', 'stim_length', 1200, 'd_save', '/Volumes/Lab/Users/Nora/pillow_test');
