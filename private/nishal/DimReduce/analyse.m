startup_analyse_tenessee
addpath('/Volumes/Analysis/nishal/DataHigh1.1/');

datafile = '2012-08-09-3/data003';
datarun=load_data(datafile)
datarun=load_sta(datarun)
datarun=load_params(datarun)
datarun=load_ei(datarun,'all','array_type',519);

cell_ids = datarun.cell_ids;
numCell=length(cell_ids);
neuronFile=edu.ucsc.neurobiology.vision.io.NeuronFile('/Volumes/Analysis/2012-08-09-3/data003/data003.neurons');

vision_spk_times=cell(length(cell_ids),1);
for icell=1:length(cell_ids)
vision_spk_times{icell}= neuronFile.getSpikeTimes(cell_ids(icell));
end

numTrials=100;
[dataTrials,allElecData,trialStartIndices]=getTrialDataMultipleElectrodes('/Volumes/Data/2012-08-09-3/data003',23,11*numTrials,numTrials,[23]);

D=struct([]);
D(numTrials-1).data=[];


for iTrial=1:numTrials-1
    iTrial
data=zeros(numCell,11*20000);
for icell=1:numCell
vec=zeros(11*20000,1);
trialCellSpikeTimes= vision_spk_times{icell}(vision_spk_times{icell}>=trialStartIndices(iTrial) & vision_spk_times{icell}<trialStartIndices(iTrial+1)) - trialStartIndices(iTrial);
trialCellSpikeTimes = trialCellSpikeTimes(trialCellSpikeTimes<11*20000);
trialCellSpikeTimes=trialCellSpikeTimes(trialCellSpikeTimes>0);
vec(trialCellSpikeTimes)=1;
data(icell,:)=vec';
end
D(iTrial).data=data;
end
