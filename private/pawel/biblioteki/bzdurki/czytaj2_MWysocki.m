NS_GlobalConstants=NS_GenerateGlobalConstants(500);

%1) Reading raw data
%javaaddpath 'C:\home\pawel\praca\Vision6-std-executable\Vision.jar'; %define path to Vision jar file
full_path='I:\backup1\2012-09-27-4\data000'; %define path to raw data file
rawFile=edu.ucsc.neurobiology.vision.io.RawDataFile(full_path); 
RawData=rawFile.getData(1440000,40000)'; %the output is 513x40000 array (for 512-channel file). First index is channel number, and there is 40000 samples for each channel. The first sample is sample number 100000, as specified in the first argument.
plot(RawData(212,:)); % elektroda 211!
%break;

%2) Reading some neuron information
paramsFile=edu.ucsc.neurobiology.vision.io.ParametersFile('I:\backup1\analysis\2012-09-27-4\data000\data000.params');
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile('I:\backup1\analysis\2012-09-27-4\data000\data000.neurons');
idList = neuronFile.getIDList(); %this imports list of neurons IDs from the file that is defined above as 'neuronFile';

NeuronID=5898;
spikeTimes = neuronFile.getSpikeTimes(NeuronID)' % for given neuron, import the spikes times

%3. Plotting the primary channel
%CenterChannel=60;
N=200; % ile spikow
L=80; % ile probek - 4 ms
spikes=zeros(N,L);
for i=1:N
    t=spikeTimes(i);
    d0=rawFile.getData(t-20,L)'; % czytaj od 1 ms przed spikiem do 3 ms po spiku
    d1=d0(2:513,:); % dane dla 1 spika, wszystkie 512 elektrod    
end

% znalezc elektrode z najwiekszym spikiem - nazwa CenterChannel
% jak znalezc sasiiednie elektrody

electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(500); %define the electrode map - must be different than 1 for silicon probes
Radius=1;
ChannelsPlot=electrodeMap.getAdjacentsTo(CenterChannel,Radius)';