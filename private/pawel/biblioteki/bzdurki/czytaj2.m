NS_GlobalConstants=NS_GenerateGlobalConstants(61);

%1) Reading raw data
%javaaddpath 'C:\home\pawel\praca\Vision6-std-executable\Vision.jar'; %define path to Vision jar file
full_path='G:\analysis\slices\2013-12-15-0\data004'; %define path to raw data file
rawFile=edu.ucsc.neurobiology.vision.io.RawDataFile(full_path); 
RawData=rawFile.getData(1440000,40000)'; %the output is 65x40000 array (for 64-channel file). First index is channel number, and there is 40000 samples for each channel. The first sample is sample number 100000, as specified in the first argument.
plot(RawData(212,:));
break;

%2) Reading some neuron information
paramsFile=edu.ucsc.neurobiology.vision.io.ParametersFile('E:\analysis\2008-12-06-0\2008-12-06-0\data000\data000.params');
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile('E:\analysis\2008-12-06-0\2008-12-06-0\data000\data000.neurons');
idList = neuronFile.getIDList(); %this imports list of neurons IDs from the file that is defined above as 'neuronFile';
NeuronID=886;
spikeTimes = neuronFile.getSpikeTimes(NeuronID)'; % for given neuron, import the spikes times

%3. Plotting the primary channel
CenterChannel=60;
N=200;
L=80;
spikes=zeros(N,L);
for i=1:N
    t=spikeTimes(i);
    d0=rawFile.getData(t-20,L)';
    d1=d0(CenterChannel+1,:); % add 1, since the first channel in the raw data is the TTL channel and not electrode 1
    spikes(i,:)=d1;
end
size(d0)
size(d1)
figure(1)
plot(spikes')
grid on;

%4. More electrodes

electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(1); %define the electrode map - must be different than 1 for silicon probes
Radius=1;
ChannelsPlot=electrodeMap.getAdjacentsTo(CenterChannel,Radius)';

spikes2=zeros(N,numel(ChannelsPlot),L); %spike,channel,sample
for i=1:N
    t=spikeTimes(i);
    d0=rawFile.getData(t-20,L)';
    d1=d0(ChannelsPlot+1,:);
    spikes2(i,:,:)=d1;
end
FigureProperties=struct('FigureNumber',2,'Subplot',[2 3 3],'TimeRange',[1 55],'AmplitudeRange',[-400 -300],'FontSize',24,'Colors',['k' 'r' 'b' 'k' 'g' 'm' 'c'],'LineWidth',1,'YLabel','input signal [\muV]');
y=NS_PlotClustersOfSignaturesOnArrayLayout(spikes2,ChannelsPlot,ones(1,N),1,FigureProperties,NS_GlobalConstants);

%5. Calculate EI
EI=mean(spikes2); %that was easy

%6. 
FigureProperties=struct('FigureNumber',3,'Subplot',[2 3 3],'TimeRange',[1 55],'AmplitudeRange',[-400 -300],'FontSize',24,'Colors',['k' 'r' 'b' 'k' 'g' 'm' 'c'],'LineWidth',2,'YLabel','input signal [\muV]');
y=NS_PlotClustersOfSignaturesOnArrayLayout(EI,ChannelsPlot,[1],1,FigureProperties,NS_GlobalConstants);