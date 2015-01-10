%% add path to nora's folder for GLM code
location_of_git_repo='/home/vision/Nishal/matlab/';
addpath(genpath([location_of_git_repo '/private/nora']));

% Java library
javaaddpath('/Volumes/Lab/Development/vision7/Vision.app/Contents/Resources/Java/Vision.jar');
addpath(genpath('../../create_act2'));
%% Fit the GLM from a classification run or other white noise run
%fittedGLM=glm_fit_from_WNrun({151}, '2014-11-05-2/data009_nps', 'RGB-10-2-0.48-11111-32x32', 900, '~/Nishal/GLM_fits');
% Or load a GLM fit that already exists
 load('/Volumes/Analysis/nora/nishal_glmfits/15min/151.mat');

%% Load the rawMovie to make predictions for
testmovie_filename='/Volumes/Data/2014-11-05-2/visual/18.rawMovie';
testmovie=get_rawmovie(testmovie_filename,5760);
testmovie=permute(testmovie,[2 3 1]);

%% Load the datarun to compare (where the above movie was actually run)
datarun=load_data('2014-11-05-2/data010_from_data009_nps');
datarun=load_neurons(datarun);

%% Evaluate the GLM fit and plot the rasters
sim_rast_tv1=[];
sim_rast_tv2=[];

for imodel_run=1:50
    imodel_run
    % Generate response to WN 
  
    % Calculate STA 
    
    % Generate null movie from STA calculated above ? 
      null_movie_compute_ts
      testmovie_filename='~/Nishal/TS_data/18.rawMovie';
      testmovie=get_rawmovie(testmovie_filename,2880);
      testmovie=permute(testmovie,[2 3 1]);

    % Generate rasters
x=GLM_predict(fittedGLM, datarun, testmovie, 30);
%plotraster(x,fittedGLM,'labels',true,'raster_length',48)

    % Analyse Rasters.. Calculate PSTH, etc
rec_rast=x.rasters.glm_sim;
calculate_psth

cond_times3 = time>=12*2+2 & time<12*3-2;
cond_times1 = time>=12*0+2 & time<12*1-2;
 cond_times2 = time>=12*1+2 & time<12*2-2;
 cond_times4 = time>=12*3+2 & time<12*4-2;

sim_rast_tv1(imodel_run)=sqrt(var(PSTH_rec(cond_times1)));
sim_rast_tv2(imodel_run)=sqrt(var(PSTH_rec(cond_times2)));
% sim_rast_tv3(imodel_run)=sqrt(var(PSTH_rec(cond_times3)));
% sim_rast_tv4(imodel_run)=sqrt(var(PSTH_rec(cond_times4)));
end

%%

testmovie_filename='/Volumes/Data/2014-11-05-2/visual/18.rawMovie';
testmovie=get_rawmovie(testmovie_filename,5760);
testmovie=permute(testmovie,[2 3 1]);
x=GLM_predict(fittedGLM, datarun, testmovie, 30);
%plotraster(x,fittedGLM,'labels',true,'raster_length',48)

rec_rast=x.rasters.recorded;
calculate_psth
rec_rast_tv1=sqrt(var(PSTH_rec(cond_times1)));
rec_rast_tv2=sqrt(var(PSTH_rec(cond_times2)));
rec_rast_tv3=sqrt(var(PSTH_rec(cond_times3)));
rec_rast_tv4=sqrt(var(PSTH_rec(cond_times4)));

[X1,N1]=hist(sim_rast_tv1);
[X2,N2]=hist(sim_rast_tv2);
[X3,N3]=hist(sim_rast_tv3);
[X4,N4]=hist(sim_rast_tv4);

figure;
plot(N1,X1,'k');
hold on
plot([rec_rast_tv1,rec_rast_tv1],[0,10],'r');

hold on
plot(N2,X2,'b');
hold on;
plot([rec_rast_tv2,rec_rast_tv2],[0,10],'g');

hold on
plot(N3,X3,'c');
hold on;
plot([rec_rast_tv3,rec_rast_tv3],[0,10],'m');

hold on
plot(N4,X4,'y');
hold on;
plot([rec_rast_tv4,rec_rast_tv4],[0,10],'g');

legend('Sim Condition 1','Rec Condition 1','Sim Cond 2','Rec Cond 2','Sim Cond 3','Rec Cond 3','Sim Cond 4','Rec Cond 4');
