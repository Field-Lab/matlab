neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile('D:\Home\Pawel\analysis\slices\2013\2010-09-14-0\data002\data002.neurons');

NeuronID=48

SeedEl = neuronFile.getNeuronIDElectrode(NeuronID);
Electrodes=electrodeMap.getAdjacentsTo(SeedEl,1)';
    
spikeTimes = neuronFile.getSpikeTimes(NeuronID)'; % for given neuron, import the spikes times
NeuronEI=NS512_LoadEIFromRawData(full_path,spikeTimes(10:min(509,length(spikeTimes))),60,20);
    
EIforElectrodes=NeuronEI(Electrodes,:);
PP=max(EIforElectrodes')-min(EIforElectrodes');
RecordingElectrodeIndex=find(PP==max(PP));
RecordingElectrode=Electrodes(RecordingElectrodeIndex);

Pattern=339;
Movie=114;

DataPath='G:\analysis\slices\2010-09-14-0\Data_proc\data002';
ArtifactDataPath='G:\analysis\slices\2010-09-14-0\Data_proc\data002';

[DataTraces,ArtifactDataTraces,Channels]=NS_ReadPreprocessedData(DataPath,ArtifactDataPath,0,Pattern,Movie,0,0);
SDT=size(DataTraces);
DT=reshape(DataTraces(1:100,RecordingElectrode,:),SDT