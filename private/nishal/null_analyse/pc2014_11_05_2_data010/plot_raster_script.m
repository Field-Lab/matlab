

%%
neuronPairsRefVsNew = crossIdentifyNeuronIDs('/Volumes/Analysis/2014-11-05-2/data009', '/Volumes/Analysis/2014-11-05-2/data010',InterestingCell_vis_id);
ref_cells=neuronPairsRefVsNew(:,2);
%%
neuronPath = '/Volumes/Analysis/2014-11-05-2/data010/data010.neurons';
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile(neuronPath);
CellSpkTimes=neuronFile.getSpikeTimes(ref_cells(ref_cell_number));
TTL=neuronFile.getTTLTimes();

spks=CellSpkTimes;


TTLperTrial=floor(rawMovFrames/100)+1;

nTrials=length(TTL)/TTLperTrial;

spkColl=cell(nTrials,1);

condDuration=12;
samplesPerCondition=condDuration*20000;
ConditionStartTimes=[0:4]*samplesPerCondition+1;
nConditions=4;


spkCondColl=struct('spksColl',[]);
spkCondColl(nConditions).spksColl=cell(1,nTrials);

for iTrial=1:nTrials

TTLstart=TTL((iTrial-1)*TTLperTrial+1);
TTLend=TTL((iTrial)*TTLperTrial);

spkColl{iTrial}=(spks((spks>=TTLstart)&(spks<=TTLend))'-TTLstart);


for icond=1:nConditions
spkCondColl(icond).spksColl{iTrial}=spkColl{iTrial}(spkColl{iTrial}>=ConditionStartTimes(icond) & spkColl{iTrial}<ConditionStartTimes(icond+1)) - ConditionStartTimes(icond)+1;
end

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

figure;

for icond=1:3%nConditions
subplot(3,1,icond);
plotSpikeRaster(spkCondColl(icond).spksColl,'PlotType','vertline');
title(sprintf('%s: data0009 vis ID: %d Avg Spk Rate: %f',cond_str{icond},InterestingCell_vis_id(ref_cell_number),spkCondColl(icond).avgSpkRate));
end

%%
nConditions=3;
nTrials1=30;
figure;
for icond=1:nConditions
subplot(nConditions,1,icond);
[xPoints, yPoints]=plotSpikeRaster(spkCondColl(icond).spksColl,'PlotType','vertline');

plot(xPoints, yPoints+nTrials1, 'k');
spkCondColl(icond).xPoints=xPoints;
spkCondColl(icond).yPoints=yPoints;

end

col='mbkg';
figure('Color','w');
for icond=1:nConditions

xPoints = spkCondColl(icond).xPoints;
yPoints = spkCondColl(icond).yPoints;
nTrials1=max(yPoints(:));
plot(xPoints*120/20000, yPoints+(nConditions-icond)*nTrials1,col(icond));
hold on
ylim([0,nConditions*nTrials]);
title(sprintf('%s: data0009 vis ID: %d, Avg Spk rates (%0.02f,%0.02f,%0.02f) spks/sec',cond_str{icond},InterestingCell_vis_id(ref_cell_number),spkCondColl(1).avgSpkRate,spkCondColl(2).avgSpkRate,spkCondColl(3).avgSpkRate));
end