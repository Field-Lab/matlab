clear
NS_GlobalConstants=NS_GenerateGlobalConstants(512);

%1) Reading raw data
%javaaddpath 'C:\home\pawel\praca\Vision6-std-executable\Vision.jar'; %define path to Vision jar file
full_path='K:\analysis\2010-09-14-0\data000'; %define path to raw data file
rawFile=edu.ucsc.neurobiology.vision.io.RawDataFile(full_path); 
%RawData=rawFile.getData(1000000,40000)'; %the output is 65x40000 array (for 64-channel file). First index is channel number, and there is 40000 samples for each channel. The first sample is sample number 100000, as specified in the first argument.
%figure(5);
%plot(RawData(440,:),'b-');

%2) Reading some neuron information
paramsFile=edu.ucsc.neurobiology.vision.io.ParametersFile('C:\home\Pawel\nauka\analiza\SlicesTTX\2010-09-14-0\analysis_2013_07_18\VisionOutput\data000\data000.params');
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile('C:\home\Pawel\nauka\analiza\SlicesTTX\2010-09-14-0\analysis_2013_07_18\VisionOutput\data000\data000.neurons');
idList = neuronFile.getIDList(); %this imports list of neurons IDs from the file that is defined above as 'neuronFile';
NeuronID=6998;
spikeTimes = neuronFile.getSpikeTimes(NeuronID)'; % for given neuron, import the spikes times

%3. Plotting the primary channel
CenterChannel=336;
N=length(spikeTimes);
L=600;
spikes=zeros(N,L);
for i=1:0
    i
    t=spikeTimes(i);
    d0=rawFile.getData(t-200,L)';
    d1=d0(CenterChannel+1,:); % add 1, since the first channel in the raw data is the TTL channel and not electrode 1
    spikes(i,:)=d1;
end
size(d0)
size(d1)
figure(11)
plot(spikes')
grid on;
%break
%4. More electrodes

electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(500); %define the electrode map - must be different than 1 for silicon probes
Radius=1;
ChannelsPlot=[1:512]%electrodeMap.getAdjacentsTo(CenterChannel,Radius)';

Patterns=[339 371 409 395 463 491 32 64 47 30 45 356];
ChannelsPlot=Patterns;

spikes2=zeros(N,numel(ChannelsPlot),L); %spike,channel,sample
for i=1:N
    i
    t=spikeTimes(i);
    d0=rawFile.getData(t-100,L)';
    d1=d0(ChannelsPlot+1,:);
    spikes2(i,:,:)=d1;
end
FigureProperties=struct('FigureNumber',2,'Subplot',[2 3 3],'TimeRange',[1 55],'AmplitudeRange',[-450 -250],'FontSize',24,'Colors',['k' 'r' 'b' 'k' 'g' 'm' 'c'],'LineWidth',1,'YLabel','input signal [\muV]');
%y=NS_PlotClustersOfSignaturesOnArrayLayout(spikes2,ChannelsPlot,ones(1,N),500,FigureProperties,NS_GlobalConstants);

%5. Calculate EI
EI=mean(spikes2); %that was easy

EI2=reshape(EI,12,600);
for i=1:12
    subplot(3,4,i);
    plot(EI2(i,:));
end

break;
%6. 
FigureProperties=struct('FigureNumber',3,'Subplot',[2 3 3],'TimeRange',[1 400],'AmplitudeRange',[-340 -310],'FontSize',16,'Colors',['k' 'r' 'b' 'k' 'g' 'm' 'c'],'LineWidth',2,'YLabel','input signal [\muV]');
y=NS_PlotClustersOfSignaturesOnArrayLayout(EI,ChannelsPlot,[1],500,FigureProperties,NS_GlobalConstants);