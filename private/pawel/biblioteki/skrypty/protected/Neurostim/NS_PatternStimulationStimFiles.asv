NumberOfSamples=30000;

% * * * * * 1. Define electrodes stimulation thresholds

%NeuronIDs = [1 2 3 4 5 6 7 8 9 10];
StimElectrodes = [29 28 19 6 43 59 52 61 46 12];
StimAmps = [1 1 2.5 1.81 1.20 1 1.30 1.81 2.0 1]; % in uA, peak cathodic
NumberOfElectrodes=length(StimElectrodes);
electrodes=1:64;
Array=zeros(64);
for i = 1:NumberOfElectrodes
    Array(StimElectrodes(i),StimElectrodes(i)) = StimAmps(i);
end

% * * * * * 2. Define the pattern - choose a, b or c
%a) arbitrary waveform
Times=[1200 3400 5600 6000 18000 25000];
Patterns=StimElectrodes([1 4 5 4 2 1]);

%b) read sparse array
ArrayFileName='C:\pawel\pliki\nauka\matlab\TIMERASTER2006010601data00060mins133Neurons1p2ms';
a=impordtdata(ArrayFileName);

Pattern=[];
Times=[];
for i=1:length(StimElectrodes)
    times=find(a(1,:)~=0);
    Times=[Times times];
    patterns=ones(1,length(times))*StimElectrodes(i);
    Pattern=[Patterns patterns];
end

%c) read Vision files
paramsFile=edu.ucsc.neurobiology.vision.io.ParametersFile('C:\praca\analiza\2008-02-04-0\outputs\data005\data005000\data005000.params');
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile('C:\praca\analiza\2008-02-04-0\outputs\data005\data005000\data005000.neurons');
idList = neuronFile.getIDList();
NeuronID=886;
spikeTimes = neuronFile.getSpikeTimes(NeuronID)';



% * * * * * 3. Generate files
MovieChunk=NS_MovieChunkGenerationForExperiment(Times,NumberOfSamples,Patterns);
MovieChunksFile=[1 MovieChunk];

cd C:\pawel\pliki\nauka\matlab\; 
fid = fopen('PatternStim_electrodes','wb')
fwrite(fid,electrodes,'int32');
fclose(fid);

fid = fopen('PatternStim_patterns','wb','ieee-le.l64')
fwrite(fid,Array,'double');
fclose(fid);

fid = fopen('PatternStim_movie','wb','ieee-le.l64')
fwrite(fid,MovieChunksFile,'int32');
fclose(fid); 