
function [spkColl,spkCondColl]=plot_raster_script_withGLM_pred(datarun,WN_datafile,WN_datafile_full,Null_datafile,InterestingCell_vis_id,imov,ref_cell_number,nConditions,condDuration,cond_str,PlotConds)
%% Map neurons. Might be bad! Do better! Replace with CBP spike sorting.
neuronPairsRefVsNew = crossIdentifyNeuronIDs(WN_datafile_full, Null_datafile,InterestingCell_vis_id);
ref_cells=neuronPairsRefVsNew(:,2);
%ref_cells=InterestingCell_vis_id;
%% Load spikes
neuronPath = [Null_datafile,sprintf('/data010.neurons')];
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile(neuronPath);
CellSpkTimes=neuronFile.getSpikeTimes(ref_cells(ref_cell_number));
TTL=double(neuronFile.getTTLTimes());

spks=double(CellSpkTimes);

rawMovFrames = nConditions*condDuration*120;
TTLperTrial=floor(rawMovFrames/100)+1;

nTrials=floor(length(TTL)/TTLperTrial);

spkColl=cell(nTrials,1);

%condDuration=12;
% nConditions=6;
samplesPerCondition=condDuration*20000;
ConditionStartTimes=[0:nConditions]*samplesPerCondition+1;

spkCondColl=struct('spksColl',[]);
spkCondColl(nConditions).spksColl=cell(1,nTrials);
TTLstartOld=0;
for iTrial=1:nTrials

TTLstart=TTL((iTrial-1)*TTLperTrial+1);
TTLend=TTLstart+condDuration*nConditions*20000;%TTL((iTrial)*TTLperTrial+1);
(TTLstart-TTLstartOld)/20000
spkColl{iTrial}=(spks((spks>=TTLstart)&(spks<=TTLend))'-TTLstart);


for icond=1:nConditions
    
spkCondColl(icond).spksColl{iTrial}=spkColl{iTrial}(spkColl{iTrial}>=ConditionStartTimes(icond) & spkColl{iTrial}<ConditionStartTimes(icond+1)) - ConditionStartTimes(icond)+1;
if(isempty(spkCondColl(icond).spksColl{iTrial}))
    
spkCondColl(icond).spksColl{iTrial}=1;
end
end
TTLstartOld=TTLstart;
end

% Calculate spike rate accross conditions
for icond=1:nConditions
    numSpk=0;    
    
    for iTrial=1:nTrials
    numSpk=numSpk+length(spkCondColl(icond).spksColl{iTrial});
    end
    
    numSpk=numSpk/nTrials;
    spkCondColl(icond).avgSpkRate=numSpk/condDuration;
end
% 
% figure;
% plotSpikeRaster(spkColl,'PlotType','vertline');

% figure;
% 
% for icond=1:nConditions
% subplot(nConditions,1,icond);
% plotSpikeRaster(spkCondColl(icond).spksColl,'PlotType','vertline');
% title(sprintf('%s: data012 vis ID: %d Avg Spk Rate: %f',cond_str{icond},InterestingCell_vis_id(ref_cell_number),spkCondColl(icond).avgSpkRate));
% end

%% Plot Raster
%nConditions=6;

nTrials1=nTrials;
h=figure;
for icond=1:nConditions
subplot(nConditions,1,icond);
[xPoints, yPoints]=plotSpikeRaster(spkCondColl(icond).spksColl,'PlotType','vertline');


plot(xPoints, yPoints+nTrials1, 'k');
spkCondColl(icond).xPoints=xPoints;
spkCondColl(icond).yPoints=yPoints;

end
close(gcf);

col='rkrkrk';
figure('Color','w');
for icond=1:nConditions

xPoints = spkCondColl(icond).xPoints;
yPoints = spkCondColl(icond).yPoints;
nTrials1=max(yPoints(:));
plot(xPoints*120/20000, yPoints+(nConditions-icond)*nTrials1,col(icond));
hold on
ylim([0,nConditions*nTrials]);
%title(sprintf('%s: data%03d vis ID: %d, Avg Spk rates (%0.02f,%0.02f,%0.02f %0.02f) spks/sec',cond_str{icond},imov,InterestingCell_vis_id(ref_cell_number),spkCondColl(1).avgSpkRate,spkCondColl(2).avgSpkRate,spkCondColl(4).avgSpkRate,spkCondColl(6).avgSpkRate));
end


%% 
% figure;
% plotSpikeRaster(spkColl,'PlotType','vertline');

%%
% Fit rank 2
%fittedGLM=glm_fit_from_WNrun({InterestingCell_vis_id(ref_cell_number)}, '2014-11-05-2/data009', 'RGB-10-2-0.48-11111-32x32', 900, '~/Nishal/GLM_fits');
% Load rank 1
load(sprintf('/Volumes/Analysis/nora/nishal_glmfits/30min/%d.mat',InterestingCell_vis_id(ref_cell_number)));
 
datarun=load_data('2014-11-05-2/data009');
datarun=load_neurons(datarun);

location_of_git_repo='/home/vision/Nishal/matlab/';
addpath(genpath([location_of_git_repo '/private/nora']));

testmovie_filename='/Volumes/Data/2014-11-05-2/Visual/18.rawMovie';
testmovie=get_rawmovie(testmovie_filename,5760);
testmovie=permute(testmovie,[2 3 1]);

prediction=GLM_predict(fittedGLM,testmovie, 30);

sim_Pts=cell(nConditions,1);
    block_len = size(prediction.rasters.glm_sim,2)/nConditions;
for icond=1:nConditions

[xPoints, yPoints]=plotSpikeRaster(logical(prediction.rasters.glm_sim(:,block_len*(icond-1)+1:icond*block_len)),'PlotType','vertline');
sim_Pts{icond}.xPoints=xPoints;
sim_Pts{icond}.yPoints=yPoints;

end

%h=plotraster(prediction,fittedGLM,'labels',true,'raster_length',48)
%close(h);

col='rkrkrkrk';
figure('Color','w');
ymax=0;
icolcnt=0;
for pltcnd=PlotConds
 
    h=figure('Color','w');
    icolcnt=icolcnt+1
    xpts=sim_Pts{pltcnd}.xPoints*prediction.rasters.bintime;
    validIdx=(xpts>1)&(xpts<11);
    ypts=sim_Pts{pltcnd}.yPoints;
    plot(xpts, ypts+ymax, col(icolcnt));
    xlim([1,11]);
    set(gca,'YTick',[]);
    set(gca,'YColor','w');
    set(gca,'XColor','w')
     set(gca,'XTick',[]);
     ylim([0,max(ypts)]);
     
     s=hgexport('readstyle','raster');
     hgexport(h,sprintf('/Volumes/Analysis/nishal/Presentations/Figures for EJ/Sim_cond_%d.eps',pltcnd),s);
 %    ymax=ymax+31;
       
    h=figure('Color','w');
    icolcnt=icolcnt+1
    xpts=spkCondColl(pltcnd).xPoints*1/20000;
    validIdx=(xpts>1)&(xpts<11);
    ypts=spkCondColl(pltcnd).yPoints;
    plot(xpts, ypts+ymax, col(icolcnt));
   % ymax=ymax+31;
    xlim([1,11]);
    ylim([0,max(ypts)]);
    
    set(gca,'XColor','w')
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
    set(gca,'YColor','w');
    s=hgexport('readstyle','raster');
    hgexport(h,sprintf('/Volumes/Analysis/nishal/Presentations/Figures for EJ/Rcd_cond_%d.eps',pltcnd),s);
  
end

%xlabel('Time (s)');

%ylim([0,ymax]);
end
