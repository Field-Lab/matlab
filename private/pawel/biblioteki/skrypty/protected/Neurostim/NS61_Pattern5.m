%break;
clear;
NumberOfSamples=60000;
AllElectrodes=[1:64];
Array=eye(64,64); 

Electrodes=[1 18 37 38 41 42 50 54 59];
MovieNumbers=[26 24 22 23 22 20 21 23 26]; % 15 18 25 12 12 22 13];
Amplitudes=Electrodes;

% DO NOT FORGET THIS !!!!!!!
for i=1:length(Electrodes)
    Amplitudes(i)=NS_PulseAmplitude('C:\home\pawel\2010\analysis\retina_61\2010-09-21-0\data003',Electrodes(i),Electrodes(i),MovieNumbers(i));
end

Neurons=[76 227 256 271 391 406 541 616 691 736 856 901]

if (length(Electrodes)~=length(Neurons) | length(Electrodes)~=length(Amplitudes))
    error('Lengths of Electrodes, Amplitudes and Neurons are not identical!');
end

% 1. Apply the amplitudes
for i=1:length(Electrodes)
    Array(Electrodes(i),Electrodes(i))=Amplitudes(i);    
end

figure(1);
clf;

% 2. Natural activivy under white noise
subplot(4,1,1);
%paramsFile=edu.ucsc.neurobiology.vision.io.ParametersFile('C:\home\pawel\2010\analysis\retina_61\2010-09-21-0\data002_old\data002000\data002000.params');
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile('C:\home\pawel\2010\analysis\retina_61\2010-09-21-0\data002_old\data002000\data002000.neurons');

startTime=24e4;
TimeRange=[startTime startTime+20000];
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
Chunk1=NS_MovieChunkGenerationForExperiment(80*round(Times/80),NumberOfSamples,Patterns);

% 2. Natural activity under moving bar. 
subplot(4,1,2);

%paramsFile=edu.ucsc.neurobiology.vision.io.ParametersFile('C:\home\pawel\2010\analysis\retina_61\2010-09-21-0\data002_old\data002000\data001000.params');
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile('C:\home\pawel\2010\analysis\retina_61\2010-09-21-0\data005-mapped-data002\data005-mapped-data002.neurons');

%idList = neuronFile.getIDList();

startTime=28.37e4;
TimeRange=[startTime startTime+20000];
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
Chunk2=NS_MovieChunkGenerationForExperiment(80*round(Times/80),NumberOfSamples,Patterns);

% 3. Artificial pattern 1
Times=[];
PatternsForMovie=[];
N=length(Electrodes);

figure(1);
subplot(4,1,3);

for i=1:N    
    Times0=[];
    t=-8;
    for j=1:50
        t=t+12+4*round(15*rand); %was: t=t+10+4*round(20*rand) for the 2010-08-20 experiment
        Times0=[Times0 t];       
    end
    Times1=Times0(find(Times0<999));
    length(Times1);
    Times=[Times Times1];
    PatternsForMovie=[PatternsForMovie ones(1,length(Times1))*Electrodes(i)];
end

plot(Times*20,PatternsForMovie,'bd');
NoOfPulsesForPattern1=length(PatternsForMovie)
Chunk3=NS_MovieChunkGenerationForExperiment(Times*20,NumberOfSamples,PatternsForMovie);  

% 3. Artificial pattern 2
figure(1);
subplot(4,1,4);

Times=[];
PatternsForMovie=[];
N=length(Electrodes);
for i=1:N    
    Times0=[];
    t=-4;
    for j=1:50
        t=t+8+4*round(10*rand); %was: t=t+10+4*round(20*rand) for the 2010-08-20 experiment
        Times0=[Times0 t];       
    end
    Times1=Times0(find(Times0<999));
    length(Times1);
    Times=[Times Times1];
    PatternsForMovie=[PatternsForMovie ones(1,length(Times1))*Electrodes(i)];
end

%plot(Times*20,PatternsForMovie,'bd');
%NoOfPulsesForPattern2=length(PatternsForMovie)
%Chunk4=NS_MovieChunkGenerationForExperiment(Times*20,NumberOfSamples,PatternsForMovie);  

% 4. Create files
MovieChunksFile=[3 Chunk1 Chunk2 Chunk3];
break
cd F:\StimFiles;
fid = fopen('61_Mosaic5_el','wb');
fwrite(fid,AllElectrodes,'int32');
fclose(fid);

fid = fopen('61_Mosaic5_pt','wb','ieee-le.l64');
fwrite(fid,Array,'double');
fclose(fid);

fid = fopen('61_Mosaic5_mv','wb','ieee-le.l64');
fwrite(fid,MovieChunksFile,'int32');
fclose(fid);