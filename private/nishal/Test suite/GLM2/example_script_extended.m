%% add path to nora's folder for GLM code
location_of_git_repo='/home/vision/Nishal/matlab/';
addpath(genpath([location_of_git_repo '/private/nora']));

% Java library
javaaddpath('/Volumes/Lab/Development/vision7/Vision.app/Contents/Resources/Java/Vision.jar');
addpath(genpath('../../create_act2'));
addpath(genpath('../GLM'));
%% Fit the GLM from a classification run or other white noise run
for cellID=[3152,3331,3365,3620,3637,3692,3901,3902,3903,3904,3916,4129,4246,4291,4726,4789,4921,5059,5177,5326,5581,6006,6076,6391,6541,6725,6812,6826,6829,6856,7188,7532,7533,7651,7652,7726]
try
    %fittedGLM=glm_fit_from_WNrun({cellID}, '2014-11-05-2/data009_nps', 'RGB-10-2-0.48-11111-32x32', 900, '~/Nishal/GLM_fits');
% Or load a GLM fit that already exists

 load(sprintf('/Volumes/Analysis/nora/nishal_glmfits/15min/%d.mat',cellID));

%% Load the rawMovie to make predictions for
testmovie_filename='/Volumes/Data/2014-11-05-2/Visual/18.rawMovie';
testmovie=get_rawmovie(testmovie_filename,5760);
testmovie=permute(testmovie,[2 3 1]);



%% Whole model run .. WNSTA, generate movie, and respond
% sim_rast_tv1=[];
% sim_rast_tv2=[];
% 
% for imodel_run=1:1
%     imodel_run
%     % Generate response to WN
%     WNtime=120*60*30;
%     WNmovie =double(rand(32,32,WNtime)>0.5)-0.5;
%     x=GLM_predict(fittedGLM, WNmovie, 1);
%     
%     % Calculate STA 
%     % WNmovie made 4 dimensional
%     WNmov4D = zeros(32,32,3,WNtime);
%     for iframe=1:size(WNmovie,3);
%     WNmov4D(:,:,1,iframe)=WNmovie(:,:,iframe);
%     WNmov4D(:,:,2,iframe)=WNmovie(:,:,iframe);
%     WNmov4D(:,:,3,iframe)=WNmovie(:,:,iframe);
%     end
%     mov_params.mov=WNmov4D;
%     
%     sta_params.Filtlen=30;
%     sta_params.useTrial=1;
%     
%     cell_params.binsPerFrame=10;
%     
%     response.spksGen=x.rasters.glm_sim;
%     aa=repmat([1:WNtime],[10,1]);
%     response.mov_frame_number=aa(:);
%     
%     response = calculate_sta_ts(mov_params,response,sta_params,cell_params);
%     WNSTA = response.analyse.STA;
%         %  figure
%         %  for itime=1:sta_params.Filtlen
%         %  imagesc(squeeze((WNSTA(:,:,itime)))');colormap gray
%         %  caxis([min(WNSTA(:)),max(WNSTA(:))]);
%         %  colorbar
%         %  pause(1/120)
%         %  end
% 
%     % Generate null movie from STA calculated above ? 
%       null_movie_compute_ts
%       testmovie_filename='~/Nishal/TS_data/18.rawMovie';
%       testmovie=get_rawmovie(testmovie_filename,2880);
%       testmovie=permute(testmovie,[2 3 1]);
%        
%     % Generate rasters
% x=GLM_predict(fittedGLM, testmovie, 30);
% x.rasters.recorded = x.rasters.glm_sim;
% plotraster(x,fittedGLM,'labels',true,'raster_length',24)
% 
%     % Analyse Rasters.. Calculate PSTH, etc
% rec_rast=x.rasters.glm_sim;
% calculate_psth
% 
% cond_times3 = time>=12*2+2 & time<12*3-2;
% cond_times1 = time>=12*0+2 & time<12*1-2;
%  cond_times2 = time>=12*1+2 & time<12*2-2;
%  cond_times4 = time>=12*3+2 & time<12*4-2;
% 
% sim_rast_tv1(imodel_run)=sqrt(var(PSTH_rec(cond_times1)));
% sim_rast_tv3(imodel_run)=sqrt(var(PSTH_rec(cond_times2)));
%  sim_rast_tv2(imodel_run)=sqrt(var(PSTH_rec(cond_times1))); % Dummy
%  sim_rast_tv4(imodel_run)=sqrt(var(PSTH_rec(cond_times1))); % Dummy
% end
%% Load from datarun to compare (where the above movie was actually run)
datarun=load_data('2014-11-05-2/data010_from_data009_nps');
datarun=load_neurons(datarun);

testmovie_filename='/Volumes/Data/2014-11-05-2/visual/18.rawMovie';
testmovie=get_rawmovie(testmovie_filename,5760);
testmovie=permute(testmovie,[2 3 1]);

x=GLM_predict(fittedGLM,testmovie, 30,datarun);

rec_rast=x.rasters.recorded;
calculate_psth
rec_rast_tv1=sqrt(var(PSTH_rec(cond_times1)));
rec_rast_tv2=sqrt(var(PSTH_rec(cond_times2)));
rec_rast_tv3=sqrt(var(PSTH_rec(cond_times3)));
rec_rast_tv4=sqrt(var(PSTH_rec(cond_times4)));

%% Half model run .. Movie is the one used actually in experiment.

sim_rast_tv1=[];
sim_rast_tv2=[];
sim_rast_tv3=[]; % Dummy
sim_rast_tv4=[]; % Dummy


for imodel_run=1:30
    close all
    imodel_run
testmovie_filename='/Volumes/Data/2014-11-05-2/visual/18.rawMovie';
testmovie=get_rawmovie(testmovie_filename,5760);
testmovie=permute(testmovie,[2 3 1]);
x=GLM_predict(fittedGLM, testmovie, 30,datarun);

h=plotraster(x,fittedGLM,'labels',true,'raster_length',48)

if(imodel_run==1)
title(sprintf('Cell: %d',cellID));
filename = sprintf('/Volumes/Analysis/nishal/GLM_cells/GLM_predictions/%d_raster',cellID);
print(h, '-depsc', filename);
end

% Analyse Rasters.. Calculate PSTH, etc

rec_rast=x.rasters.glm_sim;
calculate_psth

cond_times3 = time>=12*2+2 & time<12*3-2;
cond_times1 = time>=12*0+2 & time<12*1-2;
 cond_times2 = time>=12*1+2 & time<12*2-2;
 cond_times4 = time>=12*3+2 & time<12*4-2;

sim_rast_tv1(imodel_run)=sqrt(var(PSTH_rec(cond_times1)));
sim_rast_tv2(imodel_run)=sqrt(var(PSTH_rec(cond_times2)));
sim_rast_tv3(imodel_run)=sqrt(var(PSTH_rec(cond_times3)));
sim_rast_tv4(imodel_run)=sqrt(var(PSTH_rec(cond_times4)));
end

%% Compare simulated metric and actual metric ..

[X1,N1]=hist(sim_rast_tv1);
[X2,N2]=hist(sim_rast_tv2);
[X3,N3]=hist(sim_rast_tv3);
[X4,N4]=hist(sim_rast_tv4);

h=figure;
B=bar(N1,X1,'k');
ch = get(B,'child');
set(ch,'facea',.3)
hold on
plot([rec_rast_tv1,rec_rast_tv1],[0,10],'k');

hold on
B=bar(N2,X2,'b');
ch = get(B,'child');
set(ch,'facea',.3)
hold on;
plot([rec_rast_tv2,rec_rast_tv2],[0,10],'b');

hold on
B=bar(N3,X3,'r');
ch = get(B,'child');
set(ch,'facea',.3)
hold on;
plot([rec_rast_tv3,rec_rast_tv3],[0,10],'r');

hold on
B=bar(N4,X4,'g');
ch = get(B,'child');
set(ch,'facea',.3)
hold on;
plot([rec_rast_tv4,rec_rast_tv4],[0,10],'Color','g');

legend('Sim Condition 1','Rec Condition 1','Sim Cond 2','Rec Cond 2','Sim Cond 3','Rec Cond 3','Sim Cond 4','Rec Cond 4','Location','best');
title(sprintf('Cell: %d',cellID));
filename = sprintf('/Volumes/Analysis/nishal/GLM_cells/GLM_predictions/%d',cellID);
print(h, '-depsc', filename);
catch
display ('A cell did not work');    
end
end