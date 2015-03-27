% test runconemodel with a smaller matrix
cd /netapp/snle/lab/temp_matlabcodeAH/glm_AH_19/prep/conemodels
clear; close all;

%load '/netapp/snle/lab/Experiments/Array/Analysis/akheitman/NSEM_Projects/Stimuli/NSEM_eye-120-3_0-3600/blockedmovie_schemeA_8pixel.mat'
%moviematrix0 = NSEMmovie.fitmovie.movie_byblock{1}.matrix;

load test_moviematrix.mat;

binsperframe = 8;  %roughly msec bins
framedur = (1/120);
max_rstar_sec_mult = 10000/255;


% Front Pad by 30 seconds
padframes = 60;
padval = 127;
frontpad = padval*ones(size(moviematrix0,1),size(moviematrix0,2), padframes);
prepad = double(moviematrix0);
moviematrix = cat(3,frontpad, prepad);
movie_rstar_frame = framedur * max_rstar_sec_mult * moviematrix;
checkfunction.check = true;
checkfunction.plotdir = '/netapp/snle/lab/Experiments/Array/Analysis/akheitman/NSEM_Projects/Stimuli/NSEM_eye-120-3_0-3600/model1_pA_schemeA';
checkfunction.plotblockname = sprintf('fitblock%d', 1);
model_name  = 'model1';
clear moviematrix padframes padval frontpad prepad max_rstar_sec_mult moviematrix0 


pAmpmovie = runconemodel(model_name, movie_rstar_frame, framedur, binsperframe, checkfunction);






