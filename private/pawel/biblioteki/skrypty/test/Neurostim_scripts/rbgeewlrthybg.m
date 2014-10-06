clear
NS_GlobalConstants=NS_GenerateGlobalConstants(512);

%1) Reading raw data
%javaaddpath 'C:\home\pawel\praca\Vision6-std-executable\Vision.jar'; %define path to Vision jar file
full_path='G:\uncompressed\2010-09-14-0\data000'; %define path to raw data file
rawFile=edu.ucsc.neurobiology.vision.io.RawDataFile(full_path); 
RawData=rawFile.getData(300000,20000)'; %the output is 65x40000 array (for 64-channel file). First index is channel number, and there is 40000 samples for each channel. The first sample is sample number 100000, as specified in the first argument.
figure(5);
plot(RawData(503,:),'b-');
%break;

%2) Reading some neuron information
paramsFile=edu.ucsc.neurobiology.vision.io.ParametersFile('J:\sep14_spont_EIs\data000.params');
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile('J:\sep14_spont_EIs\data000.neurons');
idList = neuronFile.getIDList(); %this imports list of neurons IDs from the file that is defined above as 'neuronFile';
NeuronID=6845;
spikeTimes = neuronFile.getSpikeTimes(NeuronID)'; % for given neuron, import the spikes times

%3. Plotting the primary channel
CenterChannel=372;
N=100;
L=3000;
spikes=zeros(N,L);
for i=1:N
    t=spikeTimes(i);
    d0=rawFile.getData(t-1000,L)';
    d1=d0(CenterChannel+1,:); % add 1, since the first channel in the raw data is the TTL channel and not electrode 1
    spikes(i,:)=d1;
end
figure(11)
h=plot(spikes');
set(h,'Color','b');
grid on;
%break

g=mean(spikes,1);
figure(12)
plot(g);
