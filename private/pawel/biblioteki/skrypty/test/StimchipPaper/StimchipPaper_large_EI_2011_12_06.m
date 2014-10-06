NS_GlobalConstants=NS_GenerateGlobalConstants(1);
ArrayID=1;
%ChannelsPlot=[8 10 13 16 5 7 11 14 18 2 3 6 12 17 64 1 4 15 22 26 41 58 60 61 33 36 47 54 56 38 44 49 53 55 19 21];
ChannelsPlot=[8 10 13 16 5 7 11 14 18 2 3 6 12 17 64 1 4 58 60 61 54 56 53 55 19 49 47 41 15 44];
ElectrodeMap([1:6 16 7 8 10:11 18 12:15 17 27 19:24 26 28 29:36 37 38:44 45 46:50 51:54 55:56 60 58:59 61:64])=[1:61]; %this is the same map as in the function NS_retina61_stim_eff_map_2011_12_02. This is done to synchronize electrode IDs in the two plotS - EI and thresholds plot.

%1) Reading raw data
full_path='K:\2010-09-21-0\data002\data002000.bin'; %define path to raw data file
rawFile=edu.ucsc.neurobiology.vision.io.RawDataFile(full_path); 

%2) Reading some neuron information
paramsFile=edu.ucsc.neurobiology.vision.io.ParametersFile('K:\Stimchip_Paper_data\data002_old\data002000\data002000.params');
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile('K:\Stimchip_Paper_data\data002_old\data002000\data002000.neurons');
idList = neuronFile.getIDList(); %this imports list of neurons IDs from the file that is defined above as 'neuronFile';
NeuronID=227;
spikeTimes = neuronFile.getSpikeTimes(NeuronID)'; % for given neuron, import the spikes times

electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(ArrayID); %define the electrode map - must be different than 1 for silicon probes

N=1000;
L=80;
spikes2=zeros(N,numel(ChannelsPlot),L); %spike,channel,sample
for i=1:N
    t=spikeTimes(i);
    d0=rawFile.getData(t-23,L)';
    d1=d0(ChannelsPlot+1,:);
    spikes2(i,:,:)=d1/0.27;
end
FigureProperties=struct('FigureNumber',2,'Subplot',[2 3 3],'TimeRange',[1 55],'AmplitudeRange',[-500 -300],'FontSize',24,'Colors',['k' 'r' 'b' 'k' 'g' 'm' 'c'],'LineWidth',1,'YLabel','input signal [\muV]');
%y=NS512_PlotClustersOfSignaturesOnArrayLayoutWithMarks(spikes2,ChannelsPlot,ones(1,N),ArrayID,FigureProperties,NS_GlobalConstants,[307],[307]);

%5. Calculate and plot EI
EI=mean(spikes2); %that was easy
for i=1:length(ChannelsPlot)
    EI(1,i,:)=EI(1,i,:)-EI(1,i,1);
end
FigureProperties=struct('FigureNumber',3,'Subplot',[2 3 3],'TimeRange',[0 60],'AmplitudeRange',[-300 150],'FontSize',20,'Colors',['k' 'r' 'r' 'k' 'g' 'm' 'c'],'LineWidth',2,'YLabel','Signal [\muV]');
t=[-0.5:0.05:3.45];
LabeledChannels=[];
%LabeledChannels=[8 10 13 16 5 7 11 14 18 2 3 6 12 17 63 64 1 4 15 22 26 41 58 60 61 62 33 36 47 54 56 59 38 44 49 53 55 19 21 28 35];

y=NS512_PlotEIOnArrayLayoutStimchipPaper2011_12_06(EI,ChannelsPlot,[1],ArrayID,FigureProperties,NS_GlobalConstants,[16 13 14 1 60 56 55],ElectrodeMap,t);
%break;
h=gcf;
FullName=['C:\home\Pawel\nauka\StimchipPaper\2012grudzien\obrazki\robocze\EI_large'];
set(h,'PaperUnits','inches');
set(h,'PaperSize',[11 12]);
set(h,'PaperPosition',[0 0 11 12]);  
print(h, '-dtiff', '-r400', FullName);