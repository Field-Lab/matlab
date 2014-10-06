NS_GlobalConstants=NS_GenerateGlobalConstants(61);

%1. Spontaneous data
full_path='J:\2008-12-06-0\data000\data000000.bin'; %define path to raw data file
rawFile=edu.ucsc.neurobiology.vision.io.RawDataFile(full_path); 
paramsFile=edu.ucsc.neurobiology.vision.io.ParametersFile('H:\analysis\2008-12-06-0\data000\data000000\data000000.params');
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile('H:\analysis\2008-12-06-0\data000\data000000\data000000.neurons');
idList = neuronFile.getIDList(); %this imports list of neurons IDs from the file that is defined above as 'neuronFile';
NeuronID=886;
spikeTimes = neuronFile.getSpikeTimes(NeuronID)'; % for given neuron, import the spikes times

CenterChannel=16;
N=50;
L=70;
spikes=zeros(N,L);
for i=1:N
    t=spikeTimes(i);
    d0=rawFile.getData(t,L)';
    d1=d0(CenterChannel+1,:); % add 1, since the first channel in the raw data is the TTL channel and not electrode 1
    spikes(i,:)=d1;
end
s1=mean(spikes); %EI of spontaneous activity
s1=s1-mean(s1);

ChannelsPlot=[1 64 63 58 60 61 62 47 54 56 59 44 49 53 55 43 46 50 52];
spikes2=zeros(N,numel(ChannelsPlot),L); %spike,channel,sample
for i=[27:100]
    t=spikeTimes(i);
    d0=rawFile.getData(t-11,L)';
    d1=d0(ChannelsPlot+1,:);
    spikes2(i,:,:)=d1;
end
EI0=mean(spikes2);
sEI0=size(EI0);
for i=1:sEI0(2)
    EI0(1,i,:)=EI0(1,i,:)-mean(EI0(1,i,:))
end

s(2,:,:)=EI0;

%2. Stimulation data
[EI,Channels]=NS512_EI_FromClusteredData('H:\analysis\2008-12-06-0\data005','H:\analysis\2008-12-06-0\data005',0,'H:\analysis\2008-12-06-0\data005\ClusterFile_005cp4',54,43,100,0);
sEI=size(EI);
for i=1:sEI(2)
    EI(1,i,:)=EI(1,i,:)-mean(EI(1,i,:))
end

s(1,:,:)=EI(1,ChannelsPlot,:);
 
FigureProperties=struct('FigureNumber',2,'Subplot',[2 3 3],'TimeRange',[1 55],'AmplitudeRange',[-50 20],'FontSize',24,'Colors',['k' 'r' 'b' 'k' 'g' 'm' 'c'],'LineWidth',2,'YLabel','signal [\muV]');
y=NS_PlotClustersOfSignaturesOnArrayLayout(s,ChannelsPlot,[1 2],1,FigureProperties,NS_GlobalConstants);