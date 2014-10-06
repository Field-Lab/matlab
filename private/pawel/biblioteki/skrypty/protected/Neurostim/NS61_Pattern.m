clear;
NumberOfSamples=60000;
AllElectrodes=[1:64];
Array=eye(64,64); 

Electrodes=[5 10 17 19 27];
Amplitudes=[2.0 1.5 1.1 2.9 1.1];
Neurons=[61 136 241 271 391];

if (length(Electrodes)~=length(Neurons) | length(Electrodes)~=length(Amplitudes))
    error('Lengths of Electrodes, Amplitudes and Neurons are not identical!');
end

% 1. Apply the amplitudes
for i=1:length(Electrodes)
    Array(Electrodes(i),Electrodes(i))=Amplitudes(i);    
end

% 2. Natural activity under moving bar. Chunk 1 - real spike times, Chunk2
% - spike times rouinded with resolution of 4 ms
figure(1);
clf;
subplot(1,2,1);

paramsFile=edu.ucsc.neurobiology.vision.io.ParametersFile('C:\home\pawel\2010\analysis\retina_61\2010-08-20-0\data004-mapped-data003\data004-mapped-data003.params');
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile('C:\home\pawel\2010\analysis\retina_61\2010-08-20-0\data004-mapped-data003\data004-mapped-data003.neurons');

idList = neuronFile.getIDList();
startTime=13.82e5;
TimeRange=[startTime startTime+20000];
%TimeRange=[1e7 4e7];

Times=[];
Patterns=[];
for i=1:length(Electrodes)
    Pattern=Electrodes(i);
    spikeTimes = neuronFile.getSpikeTimes(Neurons(i))';
    spikeTimesWithinRangeIndexes=find(spikeTimes>TimeRange(1) & spikeTimes<TimeRange(2));
    SpikeTimesToApply=spikeTimes(spikeTimesWithinRangeIndexes)-TimeRange(1);
    Times=[Times SpikeTimesToApply]; 
    Patterns=[Patterns ones(1,length(SpikeTimesToApply))*Pattern];        
    plot(SpikeTimesToApply,ones(1,length(SpikeTimesToApply))*Electrodes(i),'bd');
    hold on;
end

Chunk1=NS_MovieChunkGenerationForExperiment(Times,NumberOfSamples,Patterns);
Chunk2=NS_MovieChunkGenerationForExperiment(80*round(Times/80),NumberOfSamples,Patterns);

% 3. Artificial pattern
figure(1);
subplot(1,2,2);

Times=[];
PatternsForMovie=[];
N=length(Electrodes);
for i=1:N    
    Times0=[];
    t=-10;
    for j=1:50
        t=t+12+4*round(20*rand); %was: t=t+10+4*round(20*rand) for the 2010-08-20 experiment
        Times0=[Times0 t];       
    end
    Times1=Times0(find(Times0<997));
    length(Times1);
    Times=[Times Times1];
    PatternsForMovie=[PatternsForMovie ones(1,length(Times1))*Electrodes(i)];
end
    
plot(Times,PatternsForMovie,'bd');
Chunk3=NS_MovieChunkGenerationForExperiment(Times*20,NumberOfSamples,PatternsForMovie);  

% 4. Create files
MovieChunksFile=[1 Chunk1];
%MovieChunksFile=[3 Chunk3 Chunk3b Chunk3c];

break;
cd F:\StimFiles;
fid = fopen('61_Mosaic2_el','wb');
fwrite(fid,AllElectrodes,'int32');
fclose(fid);

fid = fopen('61_Mosaic2_pt','wb','ieee-le.l64');
fwrite(fid,Array,'double');
fclose(fid);

fid = fopen('61_Mosaic2_mv','wb','ieee-le.l64');
fwrite(fid,MovieChunksFile,'int32');
fclose(fid);