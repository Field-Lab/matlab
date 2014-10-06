clear;
Electrodes=[62:8:118 66:8:122 61:8:117 65:8:121 268:8:324 264:8:320 267:8:323 263:8:319 266:8:322 262:8:318 265:8:321 261:8:317];
Electrodes=[70:8:118 66:8:114 69:8:117 65:8:113 268:8:316 272:8:320 267:8:315 271:8:319 266:8:314 270:8:318 265:8:313 269:8:317];


NS_GlobalConstants=NS_GenerateGlobalConstants(512);
ArrayID=500;
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(ArrayID);

PatternNumber=16; %el. 287
MovieNumber=125;

%PatternNumber=40; %el. 286
%MovieNumber=123;

DataPath='E:\pawel\analysis\retina\2009-11-27-0\data001\April2010\files';
ArtifactDataPath=DataPath;
[DataTraces,ArtifactDataTraces,Channels]=NS_ReadPreprocessedData(DataPath,ArtifactDataPath,0,PatternNumber,MovieNumber,100,0);

ClusterFileName='E:\pawel\analysis\retina\2009-11-27-0\data001\April2010\files\ClusterFile_001_ID_1264';
WaveformTypes=NS_ReadClusterFile(ClusterFileName,MovieNumber,PatternNumber,100);
Artifacts=find(WaveformTypes==1);
Spikes=find(WaveformTypes==2);

ArtifactTraces=DataTraces(Artifacts,:,:);
Artifact=mean(ArtifactTraces);
TracesWithoutArtifact=NS512_SubtractArtifact(DataTraces,Artifact);

ChannelsPlot=[85 1:84 86:512];
[CorrectedTraces,EI,UniSpikesIndicCorrected]=NS512_TimingsForDetectedNeuron3(TracesWithoutArtifact(Spikes,:,:),ChannelsPlot);

EI2=EI([2:85 1 86:512],:);
EI3=EI2(Electrodes,:);
WritePathFigs='C:\home\pawel\nauka\MEA2010\prezentacja\animacje\retina\stim';

OversamplinFactor=2;
M=NS_SaveMovieFromSignatureRotated(EI3,Electrodes,[1287],OversamplinFactor,ArrayID,WritePathFigs,NS_GlobalConstants);
break;

% * * * * Spontaneous EI movie * * * * 
full_path='E:\pawel\data\retina\2009-11-27-0\data000'; %define path to raw data file
rawFile=edu.ucsc.neurobiology.vision.io.RawDataFile(full_path); 
paramsFile=edu.ucsc.neurobiology.vision.io.ParametersFile('E:\pawel\analysis\retina\2009-11-27-0\data000\data000000\data000000.params');
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile('E:\pawel\analysis\retina\2009-11-27-0\data000\data000000\data000000.neurons');
idList = neuronFile.getIDList(); %this imports list of neurons IDs from the file that is defined above as 'neuronFile';
NeuronID=1264;
spikeTimes = neuronFile.getSpikeTimes(NeuronID)'; % for given neuron, import the spikes times

N=200;
L=50;
L=50;
spikes2=zeros(N,numel(ChannelsPlot),L); %spike,channel,sample
for i=1:N
    t=spikeTimes(i);
    d0=rawFile.getData(t-10,L)';
    d1=d0([1:512]+1,:);
    spikes2(i,:,:)=d1;
end

EI=mean(spikes2); %that was easy
ei1=reshape(EI,512,L);
figure(18);
plot(min(ei1'));
EI_spont=reshape(EI(1,Electrodes,:),length(Electrodes),L);

WritePathFigs='C:\home\pawel\nauka\MEA2010\prezentacja\animacje\retina\spont';
%EI_spont=rand(96,50)*100;
M=NS_SaveMovieFromSignatureRotated(EI_spont,Electrodes,[1287],OversamplinFactor,ArrayID,WritePathFigs,NS_GlobalConstants);