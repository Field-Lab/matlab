%% Relative magnitude of sub-units

%% 
location_of_git_repo='/home/vision/Nishal/matlab/';
addpath(genpath([location_of_git_repo '/private/nora']));

% Java library

addpath(('/Volumes/Lab/Users/bhaishahster/GITs/matlab/private/nishal/create_act2'));
javaaddpath('/Volumes/Lab/Development/vision7/Vision.app/Contents/Resources/Java/Vision.jar');
addpath(genpath('/Volumes/Lab/Users/bhaishahster/GITs/matlab/private/nishal/Test suite/GLM/'));
addpath(genpath('/Volumes/Lab/Users/bhaishahster/GITs/matlab/private/nishal/Test suite/GLM2/'));
addpath(genpath('/Volumes/Lab/Users/bhaishahster/GITs/matlab/private/nishal'));
addpath(genpath('/Volumes/Lab/Users/bhaishahster/GITs/matlab/code'));
addpath(genpath('/Volumes/Lab/Users/bhaishahster/GITs/matlab/private/nishal/create_act_2/'));

%% 
ifgmlm=1; % 0 if nnmf
folder = '/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data038/Off parasol/detailed/';
icellType=2;


WN_datafile = '2015-03-09-2/streamed/data038/data038';
WN_datafile_short='2015-03-09-2/streamed/data038/data038';
movie_xml = 'BW-8-2-0.48-11111-40x40';
wrong_xml = 'BW-8-2-0.48-22222-40x40';
stim_length=1800;% 

datarun=load_data(WN_datafile)
datarun=load_params(datarun)
datarun=load_sta(datarun);
datarun=load_neurons(datarun);

%% Extract Max su strength  and minimum su strength 
inSU=0;
figure;
for nSU = [2,3,4,5,6,7]
    inSU = inSU+1;
    min_log=[];
    max_log=[];
    
    for cellID = [datarun.cell_types{icellType}.cell_ids];
    
       load(strcat(folder,sprintf('CellID_%d/gmlm/Cell%d_gmlm_su_%d.mat',cellID,cellID,nSU)));
       model = fitGMLM_log{nSU};
       su_norm_log=[];
       for isu=1:nSU
       su_norm_log = [su_norm_log;norm(model.Linear.filter{isu})];
       end
       min_log = [min_log;min(su_norm_log)];
       max_log = [max_log ;max(su_norm_log)];
    end
    
   ratio_log= min_log./max_log;

subplot(2,3,inSU);
scatter(max_log,ratio_log,10);
ylim([0,1]);
title(sprintf('# SU:%d',nSU))
end

