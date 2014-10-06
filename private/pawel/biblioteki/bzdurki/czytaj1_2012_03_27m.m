NS_GlobalConstants=NS_GenerateGlobalConstants(61);

%1) Reading raw data
%javaaddpath 'C:\home\pawel\praca\Vision6-std-executable\Vision.jar'; %define path to Vision jar file
full_path='C:\home\Pawel\nauka\analiza\AndrzejKoziec\StimForVision2\'; %define path to raw data file
rawFile=edu.ucsc.neurobiology.vision.io.RawDataFile(full_path); 

%2) Reading some neuron information
paramsFile=edu.ucsc.neurobiology.vision.io.ParametersFile('C:\home\Pawel\nauka\StimchipPaper\2012grudzien\dyskusje\data002_old\data002000\data002000.params');
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile('C:\home\Pawel\nauka\StimchipPaper\2012grudzien\dyskusje\data002_old\data002000\data002000.neurons');
idList = neuronFile.getIDList(); %this imports list of neurons IDs from the file that is defined above as 'neuronFile';
NeuronID=856;
spikeTimes = neuronFile.getSpikeTimes(NeuronID)'; % for given neuron, import the spikes times

spikeTimes=[1:40:5000];
N=100;
L=40;
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(500); %define the electrode map - must be different than 1 for silicon probes
ChannelsPlot=[1:513];

spikes2=zeros(N,numel(ChannelsPlot),L); %spike,channel,sample
for i=1:N
    t=spikeTimes(i)
    d0=rawFile.getData(t,L)';
    d1=d0(ChannelsPlot,:);
    spikes2(i,:,:)=d1;
end

EI=mean(spikes2(:,:,:));
EI2=reshape(EI,513,L);
figure(1)
M1=max(EI2');
M2=min(EI2');
amp=abs(M1-M2);
plot(amp)

channel=16;
ch=spikes2(:,channel+1,:);
ch1=reshape(ch,N,L);
plot(ch1')