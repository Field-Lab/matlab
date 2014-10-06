NS_GlobalConstants=NS_GenerateGlobalConstants(500);

%1) Reading raw data
%javaaddpath 'C:\home\pawel\praca\Vision6-std-executable\Vision.jar'; %define path to Vision jar file
full_path='G:\analysis\slices\2011-06-29-1\data001'; %define path to raw data file
%full_path='J:\2011-06-29-1\data001';
%full_path='J:\2011-06-29-1\data001';

rawFile=edu.ucsc.neurobiology.vision.io.RawDataFile(full_path); 

%z artefaktami: do okolo 1,120,000
%maks =210,000,000
Length=40000;
RawData=int16(rawFile.getData(1100000,40000))'; %the output is 65x40000 array (for 64-channel file). First index is channel number, and there is 40000 samples for each channel. The first sample is sample number 100000, as specified in the first argument.
figure(1)
t=[1:Length]/20;
plot(t,RawData(353,:),'bd-');

%2) Reading some neuron information
paramsFile=edu.ucsc.neurobiology.vision.io.ParametersFile('G:\analysis\slices\2011-06-29-1\Vision_out\Vision_out.params');
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile('G:\analysis\slices\2011-06-29-1\Vision_out\Vision_out.neurons');
idList = neuronFile.getIDList(); %this imports list of neurons IDs from the file that is defined above as 'neuronFile';
NeuronID=5344;
spikeTimes = neuronFile.getSpikeTimes(NeuronID)'; % for given neuron, import the spikes times

%3. Plotting the primary channel
N=100;
L=380;
spikes=zeros(N,L);
for i=1:N
    t=spikeTimes(i+5);
    d0=rawFile.getData(t-20,L)';
    d1=d0(34,:);
    spikes(i,:)=d1;
end
figure(3)
plot(spikes')
grid on;
break
%4. More electrodes
CenterChannel=5;
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(1); %define the electrode map - must be different than 1 for silicon probes
Radius=1;
ChannelsPlot=electrodeMap.getAdjacentsTo(CenterChannel,Radius)';

spikes2=zeros(N,numel(ChannelsPlot),L); %spike,channel,sample
for i=1:N
    t=spikeTimes(i);
    d0=rawFile.getData(t-20,L)';
    d1=d0(ChannelsPlot,:);
    spikes2(i,:,:)=d1;
end
FigureProperties=struct('FigureNumber',2,'Subplot',[2 3 3],'TimeRange',[1 55],'AmplitudeRange',[-550 -250],'FontSize',24,'Colors',['k' 'r' 'b' 'k' 'g' 'm' 'c'],'LineWidth',1,'YLabel','input signal [\muV]');
y=NS_PlotClustersOfSignaturesOnArrayLayout(spikes2,ChannelsPlot,ones(1,N),1,FigureProperties,NS_GlobalConstants);

%5. Calculate EI
EI=mean(spikes2); %that was easy
%EI2=reshape(EI,numel(ChannelsPlot),l); %just because we need 2-d array

%6. 
FigureProperties=struct('FigureNumber',3,'Subplot',[2 3 3],'TimeRange',[1 55],'AmplitudeRange',[-550 -250],'FontSize',24,'Colors',['k' 'r' 'b' 'k' 'g' 'm' 'c'],'LineWidth',2,'YLabel','input signal [\muV]');
y=NS_PlotClustersOfSignaturesOnArrayLayout(EI,ChannelsPlot,[1],1,FigureProperties,NS_GlobalConstants);

%FigureProperties=struct('FigureNumber',105,'Subplot',[2 3 3],'TimeRange',[1 55],'AmplitudeRange',[-550 -250],'FontSize',24,'Colors',['k' 'r' 'b' 'k' 'g' 'm' 'c'],'LineWidth',1,'YLabel','input signal [\muV]');
%y=NS_PlotClustersOfSignaturesOnArrayLayout(spikes2,ChannelsPlot,ones(1,N),1,FigureProperties,NS_GlobalConstants);
