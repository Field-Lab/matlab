function [rec_rast_tv,sim_rast,h,x_log]=model_run_spatial_pc2015_03_09_2(cellID,GLM_fit_link,applyconeNL,spkCondColl,n_modelruns)
%% add path to nora's folder for GLM code
location_of_git_repo='/home/vision/Nishal/matlab/';
addpath(genpath([location_of_git_repo '/private/nora']));

% Java library
javaaddpath('/Volumes/Lab/Development/vision7/Vision.app/Contents/Resources/Java/Vision.jar');
addpath(genpath('../../create_act2'));
addpath(genpath('../../Test suite'));
%%

clear spatial_method
spatial_method=cell(4,1);
spatial_method{1}.struct_orig=[];
spatial_method{1}.struct_null=[];

icell_list=0;
try
    icell_list=icell_list+1
%fittedGLM=glm_fit_from_WNrun({3152,3331,3365,3620,3637,3692,3901,3902,3903,3904,3916,4129,4246,4291,4726,4789,4921,5059,5177,5326,5581,6006,6076,6391,6541,6725,6812,6826,6829,6856,7188,7532,7533,7651,7652,7726}, '2014-11-05-2/data009_nps', 'RGB-10-2-0.48-11111-32x32', 900, '/Volumes/Analysis/nora/nishal_glmfits/15min_rank2');
%save(sprintf('/Volumes/Analysis/nora/nishal_glmfits/15min_rank2/%d.mat',cellID),'fittedGLM');

load(sprintf([GLM_fit_link,'CellID_',num2str(cellID),'.mat']));

%% Replace fitted linear filter with STA - better filter?
% 
% sta_filter=fittedGLM.cellinfo.WN_STA;
% cell_params.STAlen=30;
% sta_filt{1}=sta_filter;
% [new_stas,totalMaskAccept,CellMasks]=clipSTAs(sta_filt,cell_params);
% 
% sta_filter=squeeze(sum(sta_filter,3));
% xcoords = fittedGLM.linearfilters.Stimulus.x_coord;
% ycoords = fittedGLM.linearfilters.Stimulus.y_coord;
% totalMaskAccept=totalMaskAccept(ycoords,xcoords);
% sta_filter = sta_filter(ycoords,xcoords,:).*repmat(totalMaskAccept,[1,1,30]);
% sta_filter(:,:,1:14)=0;
% fittedGLM.linearfilters.Stimulus.Filter = sta_filter*max(abs(fittedGLM.linearfilters.Stimulus.Filter(:)))/max(abs(sta_filter(:)));
h1= figure;

subplot(1,2,1);
imagesc(fittedGLM.linearfilters.Stimulus.Filter(:,:,6));
colormap gray
caxis([min(fittedGLM.linearfilters.Stimulus.Filter(:)),max(fittedGLM.linearfilters.Stimulus.Filter(:))]);
pause(1)

subplot(1,2,2);
ssta_dummy = zeros(size(fittedGLM.linearfilters.Stimulus.Filter,1)^2,30);
for itime =1:30
    xx=fittedGLM.linearfilters.Stimulus.Filter(:,:,itime);
    ssta_dummy(:,itime)=xx(:);
end

s=svd(ssta_dummy);
subplot(1,2,2);
plot(s,'*');
title(sprintf('Cell %d',cellID));


%% Test cell
WNtime=120*24;
  WNmovie =(0.48/0.5)*(double(rand(32,32,WNtime)>0.5)-0.5)+0.5;

if(applyconeNL==1)
      WNmovie=WNmovie*255;
WNmovie = movie_cone(WNmovie,1);
end

x=GLM_predict(fittedGLM, WNmovie, 50);
%plotraster(x,fittedGLM,'labels',true,'raster_length',24,'start_time',0)
figure;
plotSpikeRaster(logical(x.rasters.glm_sim),'PlotType','vertline');



 
%% Whole model run .. WNSTA, generate movie, and respond
sim_rast_tv1=[];
sim_rast_tv2=[];
x_log=cell(30,1);
for imodel_run=1:n_modelruns
    imodel_run
    % Generate response to WN
    interval=2;
    WNtime=120*60*30/interval;
    WNmovie =(0.48/0.5)*(double(rand(32,32,WNtime)>0.5)-0.5)+0.5;
    if(applyconeNL==1)
        WNmovie=WNmovie*255;
    WNmovie = movie_cone(WNmovie,interval);
    else
     WNmovie = movie_stretch (WNmovie,interval);    
    end
    
    save('~/Nishal/TS_data/WNmovie.mat','WNmovie');
    
    x=GLM_predict(fittedGLM, WNmovie, 1);
    
    % Calculate STA 
    % WNmovie made 4 dimensional
    WNmov4D = zeros(32,32,3,WNtime);
    for iframe=1:size(WNmovie,3);
    WNmov4D(:,:,1,iframe)=WNmovie(:,:,iframe);
    WNmov4D(:,:,2,iframe)=WNmovie(:,:,iframe);
    WNmov4D(:,:,3,iframe)=WNmovie(:,:,iframe);
    end
    mov_params.mov=WNmov4D;
    
    sta_params.Filtlen=30;
    sta_params.useTrial=1;
    
    cell_params.binsPerFrame=10;
    
    response.spksGen=x.rasters.glm_sim;
    aa=repmat([1:WNtime],[10,1]);
    response.mov_frame_number=aa(:);
    
    response = calculate_sta_ts(mov_params,response,sta_params,cell_params);
    WNSTA = response.analyse.STA;
%          figure
%          for itime=1:sta_params.Filtlen
%          imagesc(squeeze((WNSTA(:,:,itime)))');colormap gray
%          caxis([min(WNSTA(:)),max(WNSTA(:))]);
%          colorbar
%          pause(1/120)
%          end

    % Generate null movie from STA calculated above ? 
      mdf_file = '/Volumes/Analysis/stimuli/white-noise-xml/BW-8-2-0.48-11111-40x40.xml';
      solver_use=8;
      null_filter='STA';
      sta_spatial_method=4;
      interval=2;
      null_movie_compute_ts
      testmovie_filename='~/Nishal/TS_data/18.rawMovie';
      testmovie=get_rawmovie(testmovie_filename,raw_mov_len);
      testmovie=permute(testmovie,[2 3 1]);
      
      % Stretch frames depending on interval ..  
      if(applyconeNL==1)
          
      testmovie = movie_cone(testmovie,interval);
      else
      testmovie = movie_stretch (testmovie,interval);
      testmovie=0.48*(testmovie-mean(testmovie(:)))/(max(abs(testmovie(:)-mean(testmovie(:))))) + 0.5;
      end
      
    % Generate rasters
x=GLM_predict(fittedGLM, testmovie, length(spkCondColl(1).spksColl));
x.rasters.recorded = x.rasters.glm_sim;
%plotraster(x,fittedGLM,'labels',true,'raster_length',24)
figure;
plotSpikeRaster(logical(x.rasters.glm_sim),'PlotType','vertline');


    % Analyse Rasters.. Calculate PSTH, etc
rec_rast=x.rasters.glm_sim;
calculate_psth

cond_int = max(time)/2;
cond_times1 = time>=cond_int*0 & time<cond_int*1;
cond_times2 = time>=cond_int*1 & time<cond_int*2;

sim_rast_tv1(imodel_run)=sqrt(var(PSTH_rec(cond_times1)));
sim_rast_tv2(imodel_run)=sqrt(var(PSTH_rec(cond_times2)));

x_log{imodel_run}=x;
end

sim_rast.cond1=sim_rast_tv1;
sim_rast.cond2=sim_rast_tv2;
%% Recorded spikes and calculate PSTH
for icond =1:length(spkCondColl)
rec_rast= makeSpikeMat(spkCondColl(icond).spksColl,x.rasters.bintime,size(x.rasters.glm_sim,2)/2);
calculate_psth
rec_rast_tv(icond)=sqrt(var(PSTH_rec));
end
   
%% make figure 
close all
leg_arr=cell(length(spkCondColl)+2,1);
   col='cmkgbry';
   h=figure;
   [N1,X1] = hist(sim_rast.cond1);
   [N2,X2] = hist(sim_rast.cond2);
   plot(X1,N1,'r');
   hold on;
   leg_arr{1}='Original';
   plot(X2,N2,'k');
   leg_arr{2}='Null';
   
   hold on;
   for icond=1:length(spkCondColl)
   plot([rec_rast_tv(icond),rec_rast_tv(icond)],[0,1],col(icond));
   leg_arr{icond+2} = sprintf('Condition %d',icond);
   end
  h_leg= legend(leg_arr,'Location','best');
   set(h_leg,'FontSize',8);
catch
    display('A cell with some error');
end


end

function   testmovie_stretch = movie_stretch (testmovie,interval)
testmovie_stretch = zeros(size(testmovie,1),size(testmovie,2),size(testmovie,3)*interval);
icnt=0;
for itime=1:size(testmovie,3)
for iint = 1:interval
icnt=icnt+1;
    testmovie_stretch(:,:,icnt) = testmovie(:,:,itime);
end
end
end

function spkMat = makeSpikeMat(cell_spk,binsz,len)
nTrials=length(cell_spk);
spkMat=zeros(nTrials,len);
bintime=[0:len]*binsz;
for itrial=1:nTrials

    for itime=1:len
        if(sum(cell_spk{itrial}>=bintime(itime)*20000 & cell_spk{itrial}<bintime(itime+1)*20000)>0)
            spkMat(itrial,itime)=1;
        end
    end
end

end