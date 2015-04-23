addpath(genpath('../null_analyse/'));
addpath(genpath('../null_analyse/analyse_functions'));
%startup_null_analyse_tenessee
startup_null_analyse_bertha

%%
% Condition strings
nConditions=6;
condDuration=1272/6;
cond_str=cell(6,1);
cond_str{1}='Original';
cond_str{2}='On selected ';
cond_str{3}='All ON';
cond_str{4}='Off selected ';
cond_str{5}='All OFF';
cond_str{6}='SBC';

interestingConditions=[1,2,3,4,5,6];
%% Load Movies

rawMovFrames=1272/(2);
figure;
icnt=0;
% make pixel histogram
for imov=[1,2,4,6,8,10]
    [stim,height,width,header_size] = get_raw_movie(sprintf('/Volumes/Data/2015-03-09-2/Visual/pc2015_03_09_2_data038/%d.rawMovie',imov),rawMovFrames,1);
    subtract_movies{3}=mean(stim,1);
    subtract_movies{3}=mean(stim,1)*0+127.5;
    movie=stim-repmat(subtract_movies{3},[rawMovFrames,1,1]);
    movie=movie/255;
    
    icnt=icnt+1;
    
    subplot(3,2,icnt);
    qq=movie;
    hist(qq(:),20)
    title(sprintf('Movie pixel histogram %d',imov));
end

% make movies
interval=2;
condMov=cell(nConditions,1);
rawMovFrames=1272/(2);
icnt=0;
% make pixel histogram
for imov=[1,2,4,6,8,10]
    [stim,height,width,header_size] = get_raw_movie(sprintf('/Volumes/Data/2015-03-09-2/Visual/pc2015_03_09_2_data038/%d.rawMovie',imov),rawMovFrames,1);
    subtract_movies{3}=mean(stim,1);
    subtract_movies{3}=mean(stim,1)*0+127.5;
    movie=stim-repmat(subtract_movies{3},[rawMovFrames,1,1]);
    movie=movie/255;
    
    icnt=icnt+1;
    qq=permute(movie,[2,3,1]);
    ifcnt = 0;
    condMov{icnt}=zeros(size(qq,1),size(qq,2),size(qq,3)*interval);
    for iframe=1:size(qq,3)
        for irepeat=1:interval
            ifcnt=ifcnt+1;
            condMov{icnt}(:,:,ifcnt)=qq(:,:,iframe)+0.5; % cond mov is between 0 and 1 now!
        end
        
    end
    
end

% make contrast map
rawMovFrames=1272/(2);
figure;
icnt=0;
cMap = cell(6,1);
h=figure('Color','w');
for imov=[1,2,4,6,8,10]
    [stim,height,width,header_size] = get_raw_movie(sprintf('/Volumes/Data/2015-03-09-2/Visual/pc2015_03_09_2_data038/%d.rawMovie',imov),rawMovFrames,1);
    subtract_movies{3}=mean(stim,1);
    subtract_movies{3}=mean(stim,1)*0+127.5;
    movie=stim-repmat(subtract_movies{3},[rawMovFrames,1,1]);
    movie=movie/255;
    
    icnt=icnt+1;
    
    subplot(3,2,icnt);
    qq=movie;
    
    cMap{icnt}=contrastMap(qq);
    
    imagesc(cMap{icnt});
    caxis([3,6]);
    colorbar
    axis image
    title(sprintf('cMap: %d',imov));
end

s=hgexport('readstyle','cMap');
hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_2/data038/cMap.eps'),s);


% make cone movies
interval=2;
rawMovFrames=1272/(interval);
icnt=0;
coneMov=cell(1,1);
h=figure('Color','w');
for imov=[1,2,4,6,8,10]
    [stim,height,width,header_size] = get_raw_movie(sprintf('/Volumes/Data/2015-03-09-2/Visual/pc2015_03_09_2_data038/%d.rawMovie',imov),rawMovFrames,1);
    subtract_movies{3}=mean(stim,1);
    subtract_movies{3}=mean(stim,1)*0+127.5;
    movie=stim-repmat(subtract_movies{3},[rawMovFrames,1,1]);
    movie=movie/255;
    
    icnt=icnt+1;
    
    subplot(3,2,icnt);
    qq=permute(movie,[2,3,1]);
    
    coneMov{icnt} = movie_cone(qq,interval);
% 
%     figure;
%     for itime=1:size(qq,3)
%     subplot(1,2,1);
%     imagesc(qq(:,:,itime));
%     colormap gray
%     colorbar
%     caxis([-0.5,0.5]);
%     
%     subplot(1,2,2);
%     imagesc(coneMov{icnt}(:,:,itime));
%     colormap gray
%     colorbar
%     caxis([min(coneMov{icnt}(:)) , max(coneMov{icnt}(:))]);
%     pause(1/120);
%     end
end

%% data041 from data031


WN_datafile = '2015-03-09-2/data031/data031';
Null_datafile = '/Volumes/Analysis/2015-03-09-2/data041-from-data031_nps';
%WN_mov='/Volumes/Analysis/stimuli/white-noise-xml/RGB-10-2-0.48-11111.xml';
WN_datafile_full = '/Volumes/Analysis/2015-03-09-2/data031/data031';
neuronPath = [Null_datafile,sprintf('/data041-from-data031_nps.neurons')];
imov=14;

datarun=load_data(WN_datafile)
datarun=load_params(datarun)


cellTypeId=[1]; % 1 for On Parasols, 2 for Off parasols
InterestingCell_vis_id=[];
for icellType=cellTypeId
    icellType
    InterestingCell_vis_id=[InterestingCell_vis_id,datarun.cell_types{icellType}.cell_ids];
end

cellTypeUsed=cellTypeId*ones(length(InterestingCell_vis_id),1);
NullCells1=[2206,4773,1008]; % cell : 6826 ?
NullCells2=datarun.cell_types{1}.cell_ids;
NullCells3=[2461,1021,7188,767,4502,5911];
NullCells4=datarun.cell_types{2}.cell_ids;
NullCells5=[5732];

condDuration=10.6;
nConditions=6;
for ref_cell_number=1:length(InterestingCell_vis_id); %11
    close all
    cellID=InterestingCell_vis_id(ref_cell_number);
    % [spkColl,spkCondColl,h]=plot_raster_script_pc2015_02_24_2(datarun,WN_datafile,WN_datafile_full,Null_datafile,InterestingCell_vis_id,imov,ref_cell_number,nConditions,condDuration,cond_str,neuronPath);
    [spkColl,spkCondColl,h]=plot_raster_script_pc2015_03_09_2_light(cellID,nConditions,condDuration,cond_str,neuronPath);
   
    % plot_mosaic_pc2015_03_09_2(datarun,InterestingCell_vis_id,ref_cell_number,NullCells1,NullCells2);
   
    InterestingCell_vis_id(ref_cell_number)
    
    if(~isdir(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_2/data041/CellType_%s/CellID_%d',datarun.cell_types{cellTypeId}.name,InterestingCell_vis_id(ref_cell_number))))
        mkdir(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_2/data041/CellType_%s/CellID_%d',datarun.cell_types{cellTypeId}.name,InterestingCell_vis_id(ref_cell_number)));
    end
   % s=hgexport('readstyle','ras_mos4');
   % hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_2/data041/CellType_%s/CellID_%d.eps',datarun.cell_types{cellTypeId}.name,InterestingCell_vis_id(ref_cell_number)),s);
       
   event_tol=0.09;
   thr=0.4;
    [spkCondColl,eventOverlap,h_event] = event_count(spkCondColl,condDuration,thr,event_tol); 
   s=hgexport('readstyle','event');
   hgexport(h_event,sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_2/data041/CellType_%s/CellID_%d/events_thr_%0.02f_event_tol_%0.02f.eps',datarun.cell_types{cellTypeId}.name,InterestingCell_vis_id(ref_cell_number),thr,event_tol),s);
   save(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_2/data041/CellType_%s/CellID_%d/events_thr_%0.02f_event_tol_%0.02f.mat',datarun.cell_types{cellTypeId}.name,InterestingCell_vis_id(ref_cell_number),thr,event_tol),'eventOverlap')
    
    event_tol=0.15;
   thr=0.4;
    [spkCondColl,eventOverlap,h_event] = event_count(spkCondColl,condDuration,thr,event_tol); 
   s=hgexport('readstyle','event');
   hgexport(h_event,sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_2/data041/CellType_%s/CellID_%d/events_thr_%0.02f_event_tol_%0.02f.eps',datarun.cell_types{cellTypeId}.name,InterestingCell_vis_id(ref_cell_number),thr,event_tol),s);
   save(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_2/data041/CellType_%s/CellID_%d/events_thr_%0.02f_event_tol_%0.02f.mat',datarun.cell_types{cellTypeId}.name,InterestingCell_vis_id(ref_cell_number),thr,event_tol),'eventOverlap')

   
   %testsuite_prediction
    %[timeLogData,psthData] = psth_variability(spkCondColl,nConditions,condMovies,cond_str,InterestingCell_vis_id,imov,ref_cell_number,interestingConditions);
    % pause
end


%% data042 from data038


WN_datafile = '2015-03-09-2/streamed/data038/data038';
Null_datafile = '/Volumes/Analysis/2015-03-09-2/data042-from-data038_streamed_nps';
%WN_mov='/Volumes/Analysis/stimuli/white-noise-xml/RGB-10-2-0.48-11111.xml';
WN_datafile_full = '/Volumes/Analysis/2015-03-09-2/data038/data038';
neuronPath = [Null_datafile,sprintf('/data042-from-data038_streamed_nps.neurons')];


datarun=load_data(WN_datafile)
datarun=load_params(datarun)


cellTypeId=[1]; % 1 for On Parasols, 2 for Off parasols
InterestingCell_vis_id=[];
for icellType=cellTypeId
    icellType
    InterestingCell_vis_id=[InterestingCell_vis_id,datarun.cell_types{icellType}.cell_ids];
end

cellTypeUsed=cellTypeId*ones(length(InterestingCell_vis_id),1);
NullCells1=[2206,4548,6826,662]; % cell : 6826 ?
NullCells2=datarun.cell_types{1}.cell_ids;
NullCells3=[2596,1232,7172,842,4547,5911];
NullCells4=datarun.cell_types{2}.cell_ids;
NullCells5=[5768];

% InterestingCell_vis_id = [4548]%[662,6826,4548] %NullCells1;
condDuration=10.6;
nConditions=6;
GLM_fit_link= '/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_2/data038/CellType_ON parasol/';
for ref_cell_number=1:length(InterestingCell_vis_id); %11
    close all
    
    
    cellID=InterestingCell_vis_id(ref_cell_number);
    
    % make directory
       if(~isdir(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_2/data042/CellType_%s/CellID_%d/',datarun.cell_types{cellTypeId}.name,InterestingCell_vis_id(ref_cell_number))))
        mkdir(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_2/data042/CellType_%s/CellID_%d/',datarun.cell_types{cellTypeId}.name,InterestingCell_vis_id(ref_cell_number)));
       end
       
    % Plot recorded raster
    % [spkColl,spkCondColl,h]=plot_raster_script_pc2015_02_24_2(datarun,WN_datafile,WN_datafile_full,Null_datafile,InterestingCell_vis_id,imov,ref_cell_number,nConditions,condDuration,cond_str,neuronPath);
    [spkColl,spkCondColl,h]=plot_raster_script_pc2015_03_09_2_light(cellID,nConditions,condDuration,cond_str,neuronPath);
%    plot_mosaic_pc2015_03_09_2(datarun,InterestingCell_vis_id,ref_cell_number,NullCells1,NullCells3);
%     s=hgexport('readstyle','ras_mos4');
%     hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_2/data042/CellType_%s/CellID_%d/recorded.eps',datarun.cell_types{cellTypeId}.name,InterestingCell_vis_id(ref_cell_number)),s);

   
%     % GLM predictions
%     figure;
%     [spkCondCollGLM,h2]=plot_GLM_prediction_pc2015_03_09_2(cellID,condMov,GLM_fit_link);
%     plot_mosaic_pc2015_03_09_2(datarun,InterestingCell_vis_id,ref_cell_number,NullCells1,NullCells3);
%     
%     h3=plot_record_prediction_pc2015_03_09_2(spkCondColl,spkCondCollGLM);
%     plot_mosaic_pc2015_03_09_2(datarun,InterestingCell_vis_id,ref_cell_number,NullCells1,NullCells3);
%     s=hgexport('readstyle','ras_mos4');
%     hgexport(h3,sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_2/data042/CellType_%s/CellID_%d/GLM_pred_recorded.eps',datarun.cell_types{cellTypeId}.name,InterestingCell_vis_id(ref_cell_number)),s);
   
%     % GLM prediction magnified - for contrast studies.
%       figure;
%       mov_scales=[10,50,50,10,10,1];
%     [spkCondCollGLM,h7]=plot_GLM_prediction_movie_scaled_pc2015_03_09_2(cellID,condMov,GLM_fit_link,mov_scales);
%     plot_mosaic_pc2015_03_09_2(datarun,InterestingCell_vis_id,ref_cell_number,NullCells1,NullCells3);
%     
%     h8=plot_record_prediction_pc2015_03_09_2(spkCondColl,spkCondCollGLM);
%     plot_mosaic_pc2015_03_09_2(datarun,InterestingCell_vis_id,ref_cell_number,NullCells1,NullCells3);
    
%     % Cone run
%     [rec_rast_tv,sim_rast,h,x_log]=cone_run_spatial_pc2015_03_09_2(cellID,GLM_fit_link,spkCondColl,1,WN_datafile)
%     s=hgexport('readstyle','model_run');
%     hgexport(h5,sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_2/data042/CellType_%s/CellID_%d/model_run_cone.eps',datarun.cell_types{cellTypeId}.name,InterestingCell_vis_id(ref_cell_number)),s);
%     save (sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_2/data042/CellType_%s/CellID_%d/data_cone.mat',datarun.cell_types{cellTypeId}.name,InterestingCell_vis_id(ref_cell_number)),'rec_rast_tv','sim_rast','spkCondColl','spkCondCollGLM','model_run_pred_log');
% %
   



% Model run
%     applyconeNL=0;
%     [rec_rast_tv,sim_rast,h4,model_run_pred_log]= model_run_spatial_pc2015_03_09_2(cellID,GLM_fit_link,applyconeNL,spkCondColl,30)
%     
%     s=hgexport('readstyle','model_run');
%     hgexport(h4,sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_2/data042/CellType_%s/CellID_%d/model_run.eps',datarun.cell_types{cellTypeId}.name,InterestingCell_vis_id(ref_cell_number)),s);
%     
%     save (sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_2/data042/CellType_%s/CellID_%d/data.mat',datarun.cell_types{cellTypeId}.name,InterestingCell_vis_id(ref_cell_number)),'rec_rast_tv','sim_rast','spkCondColl','spkCondCollGLM','model_run_pred_log');
%   

        
   event_tol=0.03;
   thr=0.4;
    [spkCondColl,eventOverlap,h_event] = event_count(spkCondColl,condDuration,thr,event_tol); 
   s=hgexport('readstyle','event');
   hgexport(h_event,sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_2/data042/CellType_%s/CellID_%d/events_thr_%0.02f_event_tol_%0.02f.eps',datarun.cell_types{cellTypeId}.name,InterestingCell_vis_id(ref_cell_number),thr,event_tol),s);
   save(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_2/data042/CellType_%s/CellID_%d/events_thr_%0.02f_event_tol_%0.02f.mat',datarun.cell_types{cellTypeId}.name,InterestingCell_vis_id(ref_cell_number),thr,event_tol),'eventOverlap')
    
    event_tol=0.05;
   thr=0.4;
    [spkCondColl,eventOverlap,h_event] = event_count(spkCondColl,condDuration,thr,event_tol); 
   s=hgexport('readstyle','event');
   hgexport(h_event,sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_2/data042/CellType_%s/CellID_%d/events_thr_%0.02f_event_tol_%0.02f.eps',datarun.cell_types{cellTypeId}.name,InterestingCell_vis_id(ref_cell_number),thr,event_tol),s);
   save(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_2/data042/CellType_%s/CellID_%d/events_thr_%0.02f_event_tol_%0.02f.mat',datarun.cell_types{cellTypeId}.name,InterestingCell_vis_id(ref_cell_number),thr,event_tol),'eventOverlap')

     InterestingCell_vis_id(ref_cell_number)
    %  pause
end

%% Refine data038/data042 analysis - number of trials was wrong .. so take care of that ,, 
cellID=662;
load(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_2/data042/CellType_%s/CellID_%d/data.mat','ON parasol',cellID))
nTrials=29;
sim_rast_tv1=[];
sim_rast_tv2=[];

GLM_fit_link= sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_2/data038/CellType_ON parasol/CellID_%d.mat',cellID);
load(GLM_fit_link)
for imodel_run=1:30
spks_sim = model_run_pred_log{imodel_run}.rasters.glm_sim;
spks_sim=spks_sim(1:nTrials,:);

rec_rast=spks_sim;
calculate_psth

cond_int = max(time)/2;
cond_times1 = time>=cond_int*0 & time<cond_int*1;
cond_times2 = time>=cond_int*1 & time<cond_int*2;

sim_rast_tv1(imodel_run)=sqrt(var(PSTH_rec(cond_times1)));
sim_rast_tv2(imodel_run)=sqrt(var(PSTH_rec(cond_times2)));


end

sim_rast.cond1=sim_rast_tv1;
sim_rast.cond2=sim_rast_tv2;


close all
leg_arr=cell(length(spkCondColl)+2,1);
   col='cgkmbry';
   h=figure;
   [N1,X1] = hist(sim_rast.cond1);
   [N2,X2] = hist(sim_rast.cond2);
   bar(X1,N1,'r');
   alpha(0.3)
   hold on;
   leg_arr{1}='Original';
   bar(X2,N2,'k');
   leg_arr{2}='Null';
   
   hold on;
   for icond=1:length(spkCondColl)
   plot([rec_rast_tv(icond),rec_rast_tv(icond)],[0,1],col(icond),'LineWidth',2);
   leg_arr{icond+2} = sprintf('Condition %d',icond);
   end
  h_leg= legend(leg_arr,'Location','best');
   set(h_leg,'FontSize',8);
   
       s=hgexport('readstyle','model_run');
    hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_2/data042/CellType_ON parasol/CellID_%d/model_run_29trials.eps',cellID),s);
    
    
    
%% data038,042, cone analysis 
cellID=6826;
load(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_2/data042/CellType_%s/CellID_%d/data_cone.mat','ON parasol',cellID))
nTrials=29;
sim_rast_tv1=[];
sim_rast_tv2=[];

GLM_fit_link= sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_2/data038/CellType_ON parasol/CellID_%d.mat',cellID);
load(GLM_fit_link)
for imodel_run=1:1
spks_sim = model_run_pred_log{imodel_run}.rasters.glm_sim;
spks_sim=spks_sim(1:nTrials,:);

rec_rast=spks_sim;
calculate_psth


cond_int = max(time)/2;
cond_times1 = time>=cond_int*0.1 & time<cond_int*1;
cond_times2 = time>=cond_int*1 & time<cond_int*2;

sim_rast_tv1(imodel_run)=sqrt(var(PSTH_rec(cond_times1)));
sim_rast_tv2(imodel_run)=sqrt(var(PSTH_rec(cond_times2)));


end

figure;
plotSpikeRaster(logical(spks_sim),'PlotType','vertline');

[x1,y1]=plotSpikeRaster(logical(spks_sim(:,cond_times1)),'PlotType','vertline');
[x2,y2]=plotSpikeRaster(logical(spks_sim(:,cond_times2)),'PlotType','vertline');

figure;
subplot(2,1,1);
plot(x2/1200,y2,'k');
hold on;
plot(x1/1200,y1+29,'r');
title('rasters');
ylim([0,29*2]);
xlim([0,max(x2(:))/1200]);

subplot(2,1,2);
plot(PSTH_rec(cond_times1),'r');
hold on;
plot(PSTH_rec(cond_times2),'k');
title('PSTH');

% Make histogram
   col='rbgcmyk';
   h=figure;
   plot([sim_rast_tv1(1),sim_rast_tv1(1)],[0,1],'Color',[0.2,0.5,0.7],'LineWidth',2);
   hold on;
   plot([sim_rast_tv2(1),sim_rast_tv2(1)],[0,1],'Color',[0.7,0.2,0.5],'LineWidth',2);
   leg_arr{1}='Original';
   leg_arr{2}='Null';
   hold on;
 
   for icond=1:length(spkCondColl)
     plot([rec_rast_tv(icond),rec_rast_tv(icond)],[0,1],col(icond),'LineWidth',2);
    leg_arr{icond+2} = sprintf('Condition %d',icond);
   end
    h_leg= legend(leg_arr,'Location','best');
    set(h_leg,'FontSize',8);
   
   
%% Refine data038/data042 analysis - model run to ask if linear front end was wrong. Exclude initial response!
cellID=6826;
load(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_2/data042/CellType_%s/CellID_%d/data.mat','ON parasol',cellID))
nTrials=29;
sim_rast_tv1=[];
sim_rast_tv2=[];
cutofflow=0.1;
cutoffhigh=0;
GLM_fit_link= sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_2/data038/CellType_ON parasol/CellID_%d.mat',cellID);
load(GLM_fit_link)
for imodel_run=1:30
spks_sim = model_run_pred_log{imodel_run}.rasters.glm_sim;
spks_sim=spks_sim(1:nTrials,:);

rec_rast=spks_sim;
calculate_psth
% PSTH_rec=smoothen_psth(PSTH_rec);
cond_int = max(time)/2;
cond_times1 = time>=cond_int*cutofflow & time<cond_int*(1-cutoffhigh);
cond_times2 = time>=cond_int*(1+cutofflow)& time<cond_int*(2-cutoffhigh);

sim_rast_tv1(imodel_run)=sqrt(var(PSTH_rec(cond_times1)));
sim_rast_tv2(imodel_run)=sqrt(var(PSTH_rec(cond_times2)));


end

sim_rast.cond1=sim_rast_tv1;
sim_rast.cond2=sim_rast_tv2;


close all
leg_arr=cell(length(spkCondColl)+2,1);
   col='cgkmbry';
   h=figure;
   [N1,X1] = hist(sim_rast.cond1);
   [N2,X2] = hist(sim_rast.cond2);
   bar(X1,N1,'r');
   alpha(0.3)
   hold on;
   leg_arr{1}='Original';
   bar(X2,N2,'k');
   leg_arr{2}='Null';
   
   hold on;
   rec_rast_tv=zeros(length(spkCondColl),1);
   times2=time(1:length(time)/2);
   time_len = max(times2);
   
   cond_times=times2>cutofflow*time_len & times2<(1-cutoffhigh)*time_len;
   
   for icond=1:length(spkCondColl)
       rec_rast= makeSpikeMat(spkCondColl(icond).spksColl,model_run_pred_log{1}.rasters.bintime,size(model_run_pred_log{1}.rasters.glm_sim,2)/2);
    calculate_psth
   % PSTH_rec=smoothen_psth(PSTH_rec);
    rec_rast_tv(icond)=sqrt(var(PSTH_rec(cond_times)));
    plot([rec_rast_tv(icond),rec_rast_tv(icond)],[0,1],col(icond),'LineWidth',2);
    leg_arr{icond+2} = sprintf('Condition %d',icond);
   end
  h_leg= legend(leg_arr,'Location','best');
   set(h_leg,'FontSize',8);
   
       s=hgexport('readstyle','model_run');
    hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_2/data042/CellType_ON parasol/CellID_%d/model_run_29trials_cutoff_start_%f_end_%f.eps',cellID,cutofflow,cutoffhigh),s);
    
    
    %% 
    correlation=zeros(2,6,5,3);
    nullscale_list = [1,5,10,25,50,100];
    origscale_list = [1,10];
    convolve_list = [50,100,150,200,400];
    cutofflow=0.2;
    cutoffhigh=0.1;
    for iorigscale=[1]%1:2
        orig_scale=origscale_list(iorigscale);
        for inullscale=[1]%1:6
            null_scale=nullscale_list(inullscale);
            figure;
            mov_scales=[orig_scale,null_scale,null_scale,1,1,1];
            [spkCondCollGLM,h7]=plot_GLM_prediction_movie_scaled_pc2015_03_09_2(cellID,condMov,GLM_fit_link,mov_scales);
            
             h8=plot_record_prediction_pc2015_03_09_2(spkCondColl,spkCondCollGLM);
             plot_mosaic_pc2015_03_09_2(datarun,InterestingCell_vis_id,ref_cell_number,NullCells1,NullCells3);
              s=hgexport('readstyle','ras_mos4');
             hgexport(h8,sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_2/data042/CellType_%s/CellID_%d/contrast_raster_scale_orig_%d_null_%d.eps',datarun.cell_types{cellTypeId}.name,cellID,orig_scale,null_scale),s);
    
            h=  figure('Color','w')
            for icond=1:3
                for iconvolve=2%1:5
                   
                    convolve=convolve_list(iconvolve);
                    rec_rast=spkCondCollGLM{icond};
                    [PSTH_rec,time]=calculate_psth_fcn(convolve,fittedGLM,rec_rast);
                    PSTH_pred = PSTH_rec;
                    PSTH_pred=PSTH_pred/norm(PSTH_pred);
                    
                    rec_rast= makeSpikeMat(spkCondColl(icond).spksColl,model_run_pred_log{1}.rasters.bintime,size(model_run_pred_log{1}.rasters.glm_sim,2)/2);
                    [PSTH_rec,time]=calculate_psth_fcn(convolve,fittedGLM,rec_rast);
                    PSTH_rec=PSTH_rec/norm(PSTH_rec);
                     
                    times2=time(1:length(time)/2);
                    time_len = max(times2);
                    cond_times=times2>cutofflow*time_len & times2<(1-cutoffhigh)*time_len;
   
                    PSTH_pred = PSTH_pred(cond_times);
                    PSTH_rec = PSTH_rec(cond_times);
                    
%                   subplot(3,1,icond);
%                     plot(times2(cond_times),PSTH_pred,'k');
%                     hold on;
%                     plot(times2(cond_times),PSTH_rec,'r');
%                     h_legend=legend('Predicted','Recorded');
%                     set(h_legend,'FontSize',5,'Location','best');
%                     
%                     correlation(iorigscale,inullscale,iconvolve,icond) = (PSTH_pred-mean(PSTH_pred))*(PSTH_rec-mean(PSTH_rec))';
%                     title(sprintf('Scale : %d Correlation %f',mov_scales(icond),correlation(iorigscale,inullscale,iconvolve,icond)));
%                     [iorigscale,inullscale,icond,iconvolve,correlation(iorigscale,inullscale,iconvolve,icond)]
%                     

                    subplot(3,1,icond);
                    plot(xcorr((PSTH_pred-mean(PSTH_pred)),(PSTH_rec-mean(PSTH_rec))),'r');
                    h_legend=legend('Cross-Correlation');
                    set(h_legend,'FontSize',5,'Location','best');
                    
                    correlation(iorigscale,inullscale,iconvolve,icond) = (PSTH_pred-mean(PSTH_pred))*(PSTH_rec-mean(PSTH_rec))';
                    title(sprintf('Scale : %d Correlation %f',mov_scales(icond),correlation(iorigscale,inullscale,iconvolve,icond)));
                    [iorigscale,inullscale,icond,iconvolve,correlation(iorigscale,inullscale,iconvolve,icond)]
                    
                end
            end
             
    s=hgexport('readstyle','psth_conditions');
    hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_2/data042/CellType_ON parasol/CellID_%d/contrast_psth__xcorr_scale_orig_%d_null_%d',cellID,orig_scale,null_scale),s);
    
        end
    end
   
    % Experiments
%     figure;
%     plot(squeeze(correlation(2,:,2,:)))
%     legend('Condition 1','Condition 2','Condition 3');

%% Event analysis 1 - 5

% get some cell IDs
WN_datafile = '2015-03-09-2/data031/data031';
datarun=load_data(WN_datafile)
datarun=load_params(datarun)
% cellTypeId=[2]; % 1 for On Parasols, 2 for Off parasols
% InterestingCell_vis_id=[];
% for icellType=cellTypeId
%     icellType
%     InterestingCell_vis_id=[InterestingCell_vis_id,datarun.cell_types{icellType}.cell_ids];
% end
InterestingCell_vis_id = [46,106,767,1021,1066,1111,1591,1968,2461,3721,3961,4502,5296,5522,5911,6291,7188,7502,7697];

WN_datafile = '2015-03-09-2/streamed/data038/data038';
datarun=load_data(WN_datafile)
datarun=load_params(datarun)
cellTypeId=[2]; % 1 for On Parasols, 2 for Off parasols
InterestingCell_vis_id2=[];
for icellType=cellTypeId
    icellType
    InterestingCell_vis_id2=[InterestingCell_vis_id2,datarun.cell_types{icellType}.cell_ids];
end
% Map cell IDs
data1 = '/Volumes/Analysis/2015-03-09-2/data031/';
data2 = '/Volumes/Analysis/2015-03-09-2/streamed/data038/';
neuronPairsRefVsNew = crossIdentifyNeuronIDs(data1,data2, InterestingCell_vis_id,InterestingCell_vis_id2);

data1_cid = InterestingCell_vis_id;
data2_cid = neuronPairsRefVsNew(:,2);
data42_150ms=[];
data42_050ms=[];
data41_050ms=[];
data41_150ms=[];

for icell=1:length(data1_cid)
    data2_cid(icell)
x=load(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_2/data042/CellType_OFF parasol/CellID_%d/events_thr_0.40_event_tol_0.15.mat',data2_cid(icell)));
data42_150ms(icell) = x.eventOverlap(1,5);

x=load(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_2/data042/CellType_OFF parasol/CellID_%d/events_thr_0.40_event_tol_0.05.mat',data2_cid(icell)));
data42_050ms(icell) = x.eventOverlap(1,5);

x=load(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_2/data041/CellType_OFF Parasol/CellID_%d/events_thr_0.40_event_tol_0.05.mat',data1_cid(icell)));
data41_050ms(icell) = x.eventOverlap(1,5);

x=load(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_2/data041/CellType_OFF Parasol/CellID_%d/events_thr_0.40_event_tol_0.15.mat',data1_cid(icell)));
data41_150ms(icell) = x.eventOverlap(1,5);
end

figure('Color','w');
plot(data42_050ms,data41_150ms,'*');
hold on;
plot([0,1],[0,1],'g');
xlabel('NDF0');
ylabel('NDF2');


data42_030ms=[];
data41_090ms=[];
for icell=1:length(data1_cid)
    data2_cid(icell)
x=load(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_2/data042/CellType_OFF parasol/CellID_%d/events_thr_0.40_event_tol_0.03.mat',data2_cid(icell)));
data42_030ms(icell) = x.eventOverlap(1,5);

x=load(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_2/data041/CellType_OFF Parasol/CellID_%d/events_thr_0.40_event_tol_0.09.mat',data1_cid(icell)));
data41_090ms(icell) = x.eventOverlap(1,5);
end

figure('Color','w');
plot(data42_030ms,data41_090ms,'*');
hold on;
plot([0,1],[0,1],'g');
xlabel('NDF0');
ylabel('NDF2');


%% 



data1_cid = [2461,1021,7188,767,4502,5911];
data2_cid = [2596,1232,7172,842,4547,5911];
data42_150ms=[];
data42_050ms=[];
data41_050ms=[];
data41_150ms=[];

for icell=1:length(data1_cid)
    data2_cid(icell)
x=load(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_2/data042/CellType_OFF parasol/CellID_%d/events_thr_0.40_event_tol_0.15.mat',data2_cid(icell)));
data42_150ms(icell) = x.eventOverlap(1,4);

x=load(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_2/data042/CellType_OFF parasol/CellID_%d/events_thr_0.40_event_tol_0.05.mat',data2_cid(icell)));
data42_050ms(icell) = x.eventOverlap(1,4);

x=load(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_2/data041/CellType_OFF Parasol/CellID_%d/events_thr_0.40_event_tol_0.05.mat',data1_cid(icell)));
data41_050ms(icell) = x.eventOverlap(1,4);

x=load(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_2/data041/CellType_OFF Parasol/CellID_%d/events_thr_0.40_event_tol_0.15.mat',data1_cid(icell)));
data41_150ms(icell) = x.eventOverlap(1,4);
end

figure('Color','w');
plot(data42_050ms,data41_150ms,'*');
hold on;
plot([0,1],[0,1],'g');
xlabel('NDF0');
ylabel('NDF2');


data42_030ms=[];
data41_090ms=[];
for icell=1:length(data1_cid)
    data2_cid(icell)
x=load(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_2/data042/CellType_OFF parasol/CellID_%d/events_thr_0.40_event_tol_0.03.mat',data2_cid(icell)));
data42_030ms(icell) = x.eventOverlap(1,4);

x=load(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_2/data041/CellType_OFF Parasol/CellID_%d/events_thr_0.40_event_tol_0.09.mat',data1_cid(icell)));
data41_090ms(icell) = x.eventOverlap(1,4);
end

figure('Color','w');
plot(data42_030ms,data41_090ms,'*');
hold on;
plot([0,1],[0,1],'g');
xlabel('NDF0');
ylabel('NDF2');

%% Event analysis 1-3

% get some cell IDs
WN_datafile = '2015-03-09-2/data031/data031';
datarun=load_data(WN_datafile)
datarun=load_params(datarun)
% cellTypeId=[2]; % 1 for On Parasols, 2 for Off parasols
% InterestingCell_vis_id=[];
% for icellType=cellTypeId
%     icellType
%     InterestingCell_vis_id=[InterestingCell_vis_id,datarun.cell_types{icellType}.cell_ids];
% end
InterestingCell_vis_id = [303,338,961,1008,1068,1682,1757,2206,3122,3331,3796,4773,6212,6261,6661,6826,7430,7743];

WN_datafile = '2015-03-09-2/streamed/data038/data038';
datarun=load_data(WN_datafile)
datarun=load_params(datarun)
cellTypeId=[1]; % 1 for On Parasols, 2 for Off parasols
InterestingCell_vis_id2=[];
for icellType=cellTypeId
    icellType
    InterestingCell_vis_id2=[InterestingCell_vis_id2,datarun.cell_types{icellType}.cell_ids];
end
% Map cell IDs
data1 = '/Volumes/Analysis/2015-03-09-2/data031/';
data2 = '/Volumes/Analysis/2015-03-09-2/streamed/data038/';
neuronPairsRefVsNew = crossIdentifyNeuronIDs(data1,data2, InterestingCell_vis_id,InterestingCell_vis_id2);

data1_cid = InterestingCell_vis_id;
data2_cid = neuronPairsRefVsNew(:,2);
data42_150ms=[];
data42_050ms=[];
data41_050ms=[];
data41_150ms=[];

for icell=1:length(data1_cid)
    data2_cid(icell)

x=load(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_2/data042/CellType_ON parasol/CellID_%d/events_thr_0.40_event_tol_0.05.mat',data2_cid(icell)));
data42_050ms(icell) = x.eventOverlap(1,3);


x=load(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_2/data041/CellType_ON Parasol/CellID_%d/events_thr_0.40_event_tol_0.15.mat',data1_cid(icell)));
data41_150ms(icell) = x.eventOverlap(1,3);
end

figure('Color','w');
plot(data42_050ms,data41_150ms,'*');
hold on;
plot([0,1],[0,1],'g');
xlabel('NDF0');
ylabel('NDF2');


data42_030ms=[];
data41_090ms=[];
for icell=1:length(data1_cid)
    data2_cid(icell)
x=load(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_2/data042/CellType_ON parasol/CellID_%d/events_thr_0.40_event_tol_0.03.mat',data2_cid(icell)));
data42_030ms(icell) = x.eventOverlap(1,3);

x=load(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_2/data041/CellType_ON Parasol/CellID_%d/events_thr_0.40_event_tol_0.09.mat',data1_cid(icell)));
data41_090ms(icell) = x.eventOverlap(1,3);
end

figure('Color','w');
plot(data42_030ms,data41_090ms,'*');
hold on;
plot([0,1],[0,1],'g');
xlabel('NDF0');
ylabel('NDF2');


%% have spkCondColl, spike rate during spiking event of originial conditions

ref_cond=1;
null_cond=3;
null_event_rate=[];
ref_event_rate=spkCondColl(ref_cond).events_spks;
bin_select=zeros(size(spkCondColl(ref_cond).rec_rast,2),1);

for ievent=1:length(spkCondColl(ref_cond).events_list_start)
bin_start = spkCondColl(ref_cond).events_bin_start(ievent);
bin_end = spkCondColl(ref_cond).events_bin_end(ievent);
bin_select(bin_start:bin_end)=1;

null_event_rate(ievent)=sum(sum(spkCondColl(null_cond).rec_rast(:,bin_start:bin_end)))/(size(spkCondColl(null_cond).rec_rast,1)*spkCondColl(ref_cond).events_list_length(ievent));
end

figure;
hist(null_event_rate,20);
xlim([0,80])

ref_cond=3;
null_cond=1;
null_event_rate=[];
ref_event_rate=spkCondColl(ref_cond).events_spks;

for ievent=1:length(spkCondColl(ref_cond).events_list_start)
bin_start = spkCondColl(ref_cond).events_bin_start(ievent);
bin_end = spkCondColl(ref_cond).events_bin_end(ievent);

null_event_rate(ievent)=sum(sum(spkCondColl(null_cond).rec_rast(:,bin_start:bin_end)))/(size(spkCondColl(null_cond).rec_rast,1)*spkCondColl(ref_cond).events_list_length(ievent));
end
figure;
hist(null_event_rate,20);
xlim([0,80])

