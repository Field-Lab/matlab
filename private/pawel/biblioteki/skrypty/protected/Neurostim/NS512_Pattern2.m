%break;
Path_general='/Volumes/Stream-phoenix/Analysis/stim512/';

NumberOfElectrodes=512;

Path_piece='2012-09-27-4/';
Path=[Path_general Path_piece];

WhiteNoise_data='000';
MovingBar_data='001';
NaturalScence_data='002';

PathWhiteNoise=[Path 'data' WhiteNoise_data '/data' WhiteNoise_data '.neurons']
PathMovingBar=[Path 'data' MovingBar_data '-from-data' WhiteNoise_data '/data' MovingBar_data '-from-data' WhiteNoise_data '.neurons']
PathNaturalScene=[Path 'data' NaturalScence_data '-from-data' WhiteNoise_data '/data' NaturalScence_data '-from-data' WhiteNoise_data '.neurons']


NumberOfSamples=10000;
AllElectrodes=[1:512];
Array=eye(512,512); 

ThresholdFilePath1=[Path_general Path_piece 'stim_scan/thresholds_1']
fid1=fopen(ThresholdFilePath1,'r');
b=fread(fid1,inf,'double');
fclose(fid1);
NeuronInformation=reshape(b,length(b)/4,4);
Neurons1=NeuronInformation(:,1);
Electrodes1=NeuronInformation(:,3);
Amplitudes1=NeuronInformation(:,4);

ThresholdFilePath2=[Path_general Path_piece 'stim_scan/thresholds_2']
fid1=fopen(ThresholdFilePath2,'r');
b=fread(fid1,inf,'double');
fclose(fid1);
NeuronInformation=reshape(b,length(b)/4,4);
Neurons2=NeuronInformation(:,1);
Electrodes2=NeuronInformation(:,3);
Amplitudes2=NeuronInformation(:,4);

Neurons=[Neurons1' Neurons2']';
Electrodes=[Electrodes1' Electrodes2']';
Amplitudes=[Amplitudes1' Amplitudes2']';

if NumberOfElectrodes==519
    Electrodes=NS512_519InverseMess(Electrodes);
end

% 1. Apply the amplitudes
for i=1:length(Electrodes)
    Array(Electrodes(i),Electrodes(i))=Amplitudes(i);    
end

figure(101);
clf;

% 2. Natural activivy under white noise
subplot(2,1,1);
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile(PathWhiteNoise);
startTime=521846;
TimeRange=[startTime startTime+NumberOfSamples];
Times=[];
Patterns=[];
for i=1:length(Electrodes)
    Pattern=Electrodes(i);
    spikeTimes = neuronFile.getSpikeTimes(Neurons(i))';
    spikeTimesWithinRangeIndexes=find(spikeTimes>TimeRange(1) & spikeTimes<TimeRange(2));
    SpikeTimesToApply=spikeTimes(spikeTimesWithinRangeIndexes)-TimeRange(1);
    Times=[Times SpikeTimesToApply]; 
    Patterns=[Patterns ones(1,length(SpikeTimesToApply))*Pattern];        
    plot(SpikeTimesToApply,ones(1,length(SpikeTimesToApply))*i,'bd');
    hold on;
end
NoOfPulsesForWhiteNoise=length(Patterns)
%Chunk1=NS_MovieChunkGenerationForExperiment(40*round(Times/40),NumberOfSamples,Patterns);

% 2. Natural activity under moving bar. 
subplot(2,2,1);
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile(PathMovingBar);
startTime=521646-200;
TimeRange=[startTime startTime+NumberOfSamples];
Times=[];
Patterns=[];
for i=1:length(Electrodes)
    Pattern=Electrodes(i);
    spikeTimes = neuronFile.getSpikeTimes(Neurons(i))';
    spikeTimesWithinRangeIndexes=find(spikeTimes>TimeRange(1) & spikeTimes<TimeRange(2));
    SpikeTimesToApply=spikeTimes(spikeTimesWithinRangeIndexes)-TimeRange(1);
    Times=[Times SpikeTimesToApply]; 
    Patterns=[Patterns ones(1,length(SpikeTimesToApply))*Pattern];        
    plot(SpikeTimesToApply,ones(1,length(SpikeTimesToApply))*i,'bd');
    hold on;
end
NoOfPulsesForMovingBar=length(Patterns)
Chunk1=NS_MovieChunkGenerationForExperiment(40*round(Times/40),NumberOfSamples,Patterns);


subplot(2,2,2);
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile(PathMovingBar);
startTime=561618-200;
TimeRange=[startTime startTime+NumberOfSamples];
Times=[];
Patterns=[];
for i=1:length(Electrodes)
    Pattern=Electrodes(i);
    spikeTimes = neuronFile.getSpikeTimes(Neurons(i))';
    spikeTimesWithinRangeIndexes=find(spikeTimes>TimeRange(1) & spikeTimes<TimeRange(2));
    SpikeTimesToApply=spikeTimes(spikeTimesWithinRangeIndexes)-TimeRange(1);
    Times=[Times SpikeTimesToApply]; 
    Patterns=[Patterns ones(1,length(SpikeTimesToApply))*Pattern];        
    plot(SpikeTimesToApply,ones(1,length(SpikeTimesToApply))*i,'bd');
    hold on;
end
NoOfPulsesForMovingBar=length(Patterns)
Chunk2=NS_MovieChunkGenerationForExperiment(40*round(Times/40),NumberOfSamples,Patterns);


subplot(2,2,3);
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile(PathMovingBar);
startTime=201535-200;
TimeRange=[startTime startTime+NumberOfSamples];
Times=[];
Patterns=[];
for i=1:length(Electrodes)
    Pattern=Electrodes(i);
    spikeTimes = neuronFile.getSpikeTimes(Neurons(i))';
    spikeTimesWithinRangeIndexes=find(spikeTimes>TimeRange(1) & spikeTimes<TimeRange(2));
    SpikeTimesToApply=spikeTimes(spikeTimesWithinRangeIndexes)-TimeRange(1);
    Times=[Times SpikeTimesToApply]; 
    Patterns=[Patterns ones(1,length(SpikeTimesToApply))*Pattern];        
    plot(SpikeTimesToApply,ones(1,length(SpikeTimesToApply))*i,'bd');
    hold on;
end
NoOfPulsesForMovingBar=length(Patterns)
Chunk3=NS_MovieChunkGenerationForExperiment(40*round(Times/40),NumberOfSamples,Patterns);



subplot(2,2,4);
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile(PathMovingBar);
startTime=1081589-200;
TimeRange=[startTime startTime+NumberOfSamples];
Times=[];
Patterns=[];
for i=1:length(Electrodes)
    Pattern=Electrodes(i);
    spikeTimes = neuronFile.getSpikeTimes(Neurons(i))';
    spikeTimesWithinRangeIndexes=find(spikeTimes>TimeRange(1) & spikeTimes<TimeRange(2));
    SpikeTimesToApply=spikeTimes(spikeTimesWithinRangeIndexes)-TimeRange(1);
    Times=[Times SpikeTimesToApply]; 
    Patterns=[Patterns ones(1,length(SpikeTimesToApply))*Pattern];        
    plot(SpikeTimesToApply,ones(1,length(SpikeTimesToApply))*i,'bd');
    hold on;
end
NoOfPulsesForMovingBar=length(Patterns)
Chunk4=NS_MovieChunkGenerationForExperiment(40*round(Times/40),NumberOfSamples,Patterns);

% * * * * * *  Natural scene * * * * * * *

figure(200)
clf
subplot(2,2,1);
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile(PathNaturalScene);
startTime=4410194;
TimeRange=[startTime startTime+NumberOfSamples];
Times=[];
Patterns=[];
for i=1:length(Electrodes)
    Pattern=Electrodes(i);
    spikeTimes = neuronFile.getSpikeTimes(Neurons(i))';
    spikeTimesWithinRangeIndexes=find(spikeTimes>TimeRange(1) & spikeTimes<TimeRange(2));
    SpikeTimesToApply=spikeTimes(spikeTimesWithinRangeIndexes)-TimeRange(1);
    Times=[Times SpikeTimesToApply]; 
    Patterns=[Patterns ones(1,length(SpikeTimesToApply))*Pattern];        
    plot(SpikeTimesToApply,ones(1,length(SpikeTimesToApply))*i,'bd');
    hold on;
end
NoOfPulses=length(Patterns)
Chunk5=NS_MovieChunkGenerationForExperiment(40*round(Times/40),NumberOfSamples,Patterns);


subplot(2,2,2);
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile(PathNaturalScene);
startTime=1110158;
TimeRange=[startTime startTime+NumberOfSamples];
Times=[];
Patterns=[];
for i=1:length(Electrodes)
    Pattern=Electrodes(i);
    spikeTimes = neuronFile.getSpikeTimes(Neurons(i))';
    spikeTimesWithinRangeIndexes=find(spikeTimes>TimeRange(1) & spikeTimes<TimeRange(2));
    SpikeTimesToApply=spikeTimes(spikeTimesWithinRangeIndexes)-TimeRange(1);
    Times=[Times SpikeTimesToApply]; 
    Patterns=[Patterns ones(1,length(SpikeTimesToApply))*Pattern];        
    plot(SpikeTimesToApply,ones(1,length(SpikeTimesToApply))*i,'bd');
    hold on;
end
NoOfPulses=length(Patterns)
Chunk6=NS_MovieChunkGenerationForExperiment(40*round(Times/40),NumberOfSamples,Patterns);


subplot(2,2,3);
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile(PathNaturalScene);
startTime=2210225;
TimeRange=[startTime startTime+NumberOfSamples];
Times=[];
Patterns=[];
for i=1:length(Electrodes)
    Pattern=Electrodes(i);
    spikeTimes = neuronFile.getSpikeTimes(Neurons(i))';
    spikeTimesWithinRangeIndexes=find(spikeTimes>TimeRange(1) & spikeTimes<TimeRange(2));
    SpikeTimesToApply=spikeTimes(spikeTimesWithinRangeIndexes)-TimeRange(1);
    Times=[Times SpikeTimesToApply]; 
    Patterns=[Patterns ones(1,length(SpikeTimesToApply))*Pattern];        
    plot(SpikeTimesToApply,ones(1,length(SpikeTimesToApply))*i,'bd');
    hold on;
end
NoOfPulses=length(Patterns)
Chunk7=NS_MovieChunkGenerationForExperiment(40*round(Times/40),NumberOfSamples,Patterns);


subplot(2,2,4);
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile(PathNaturalScene);
startTime=3310127;
TimeRange=[startTime startTime+NumberOfSamples];
Times=[];
Patterns=[];
for i=1:length(Electrodes)
    Pattern=Electrodes(i);
    spikeTimes = neuronFile.getSpikeTimes(Neurons(i))';
    spikeTimesWithinRangeIndexes=find(spikeTimes>TimeRange(1) & spikeTimes<TimeRange(2));
    SpikeTimesToApply=spikeTimes(spikeTimesWithinRangeIndexes)-TimeRange(1);
    Times=[Times SpikeTimesToApply]; 
    Patterns=[Patterns ones(1,length(SpikeTimesToApply))*Pattern];        
    plot(SpikeTimesToApply,ones(1,length(SpikeTimesToApply))*i,'bd');
    hold on;
end
NoOfPulses=length(Patterns)
Chunk8=NS_MovieChunkGenerationForExperiment(40*round(Times/40),NumberOfSamples,Patterns);


% 4. Create files
MovieChunksFile=[8 Chunk1 Chunk2 Chunk3 Chunk4 Chunk5 Chunk6 Chunk7 Chunk8];

cd([Path 'stim_files']);

if NumberOfElectrodes==519
    ElectrodeNumberLabel='_519';
else
    ElectrodeNumberLabel='';
end

fid = fopen(['Pattern512' ElectrodeNumberLabel '_el'],'wb');
fwrite(fid,AllElectrodes,'int32');
fclose(fid);

fid = fopen(['Pattern512' ElectrodeNumberLabel '_pt'],'wb','ieee-le.l64');
fwrite(fid,Array,'double');
fclose(fid);

fid = fopen(['Pattern512' ElectrodeNumberLabel '_mv'],'wb','ieee-le.l64');
fwrite(fid,MovieChunksFile,'int32');
fclose(fid);