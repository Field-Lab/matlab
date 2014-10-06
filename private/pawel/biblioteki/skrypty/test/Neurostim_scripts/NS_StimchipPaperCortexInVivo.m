NS_GlobalConstants=NS_GenerateGlobalConstants(61);

%1. Spontaneous data
full_path='E:\pawel\data\in_vivo\2009-07-25\data003\data003000.bin'; %define path to raw data file
rawFile=edu.ucsc.neurobiology.vision.io.RawDataFile(full_path); 
paramsFile=edu.ucsc.neurobiology.vision.io.ParametersFile('E:\pawel\analysis\in-vivo\2009-07-25\data003\data003000\data003000.params');
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile('E:\pawel\analysis\in-vivo\2009-07-25\data003\data003000\data003000.neurons');
idList = neuronFile.getIDList(); %this imports list of neurons IDs from the file that is defined above as 'neuronFile';
NeuronID=228;
spikeTimes = neuronFile.getSpikeTimes(NeuronID)'; % for given neuron, import the spikes times

CenterChannel=16;
N=60;
L=80;
spikes=zeros(N,L);
for i=1:N
    t=spikeTimes(i);
    d0=rawFile.getData(t-8,L)';
    d1=d0(CenterChannel+1,:); % add 1, since the first channel in the raw data is the TTL channel and not electrode 1
    spikes(i,:)=d1;
end
s1=mean(spikes); %EI of spontaneous activity
s1=s1-mean(s1);

%2. Stimulation data
StimDataPath='E:\pawel\analysis\in-vivo\2009-07-25\data006';
[EI,Channels]=NS512_EI_FromClusteredData(StimDataPath,StimDataPath,0,'E:\pawel\analysis\in-vivo\2009-07-25\data006\ClusterFile_006',1,7,100,0);
s2=reshape(EI(1,16,:),1,80);
s2=s2-mean(s2);

%3. Plotting
%subplot('Position',[0.05 0.7 0.26 0.26]);
figure(1)
t=[1:80]/20;
h=plot(t,s1*2.2,t-0.018,s2);
set(h(1),'LineWidth',2);
set(h(2),'LineWidth',2);
axis([0 2.5 -50 20]);
grid on;
h=gca;
set(h,'FontSize',16);

%4. Stimulation data - traces plotting
[DataTraces,ArtifactDataTraces,Channels]=NS_ReadPreprocessedData(StimDataPath,StimDataPath,0,1,7,100,0);
Traces=reshape(DataTraces(:,16,:),98,80)+340;

ClusterIndex=NS_ReadClusterFile('E:\pawel\analysis\in-vivo\2009-07-25\data006\ClusterFile_006',7,1,100);

figure(2)
clf
for i=2:98
    s=Traces(i,:);
    h=plot(t,s);
    index=ClusterIndex(i)
    if index==1
        set(h,'Color','k');
    else
        set(h,'Color','r');
    end
    hold on;
end
grid on
axis([0 2.5 -100 250]);

cd E:\pawel\data\in_vivo\2009-07-25;
FileName='006';
[Pulse,Status]=NS_FindPulseShapeForMovie(FileName,16,7,NS_GlobalConstants);
%figure(3)
[Amplitude,PlotPointer]=NS_PlotStimulationPulse(Pulse,Status,16,2,'b-',NS_GlobalConstants);
set(PlotPointer,'LineWidth',2);