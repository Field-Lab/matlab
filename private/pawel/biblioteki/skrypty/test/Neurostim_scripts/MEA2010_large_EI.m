NS_GlobalConstants=NS_GenerateGlobalConstants(512);
ArrayID=500;
ChannelsPlot=[292:8:324 296:8:320 291:8:323 295:8:319 290:8:322 294:8:318 289:8:321 293:8:317];
ChannelsPlot=[300:8:316 296:8:312 299:8:315 295:8:311 298:8:314 294:8:310 297:8:313 293:8:309 289:292 317:320 285:288];


%1) Reading raw data
full_path='E:\pawel\data\retina\2009-11-27-0\data000\data000000.bin'; %define path to raw data file
rawFile=edu.ucsc.neurobiology.vision.io.RawDataFile(full_path); 

%2) Reading some neuron information
paramsFile=edu.ucsc.neurobiology.vision.io.ParametersFile('E:\pawel\analysis\retina\2009-11-27-0\data000\data000000\data000000.params');
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile('E:\pawel\analysis\retina\2009-11-27-0\data000\data000000\data000000.neurons');
idList = neuronFile.getIDList(); %this imports list of neurons IDs from the file that is defined above as 'neuronFile';
NeuronID=4546;
spikeTimes = neuronFile.getSpikeTimes(NeuronID)'; % for given neuron, import the spikes times

electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(ArrayID); %define the electrode map - must be different than 1 for silicon probes

N=500;
L=80;
spikes2=zeros(N,numel(ChannelsPlot),L); %spike,channel,sample
for i=1:N
    t=spikeTimes(i);
    d0=rawFile.getData(t-27,L)';
    d1=d0(ChannelsPlot+1,:);
    spikes2(i,:,:)=d1;
end
FigureProperties=struct('FigureNumber',2,'Subplot',[2 3 3],'TimeRange',[1 55],'AmplitudeRange',[-500 -300],'FontSize',24,'Colors',['k' 'r' 'b' 'k' 'g' 'm' 'c'],'LineWidth',1,'YLabel','input signal [\muV]');
y=NS512_PlotClustersOfSignaturesOnArrayLayoutWithMarks(spikes2,ChannelsPlot,ones(1,N),ArrayID,FigureProperties,NS_GlobalConstants,[307],[307]);

%5. Calculate and plot EI
EI=mean(spikes2); %that was easy
for i=1:length(ChannelsPlot)
    EI(1,i,:)=EI(1,i,:)-EI(1,i,1);
end
FigureProperties=struct('FigureNumber',3,'Subplot',[2 3 3],'TimeRange',[1 60],'AmplitudeRange',[-150 60],'FontSize',18,'Colors',['k' 'b' 'r' 'k' 'g' 'm' 'c'],'LineWidth',1.5,'YLabel','Signal [\muV]');
t=[-0.5:0.05:3.45];
LabeledChannels=[304 307 299 303 297 293];
y=NS512_PlotClustersOfSignaturesOnArrayLayoutWithMarksRotated(EI,ChannelsPlot,[1],ArrayID,FigureProperties,NS_GlobalConstants,[],LabeledChannels,t);

h=gcf;
FullName=['C:\home\pawel\nauka\MEA2010\prezentacja\figures_pok107\EI_large_with_labels'];
set(h,'PaperUnits','inches');
set(h,'PaperSize',[9 5]);
set(h,'PaperPosition',[0 0 9 5]);  
print(h, '-dtiff', '-r120', FullName);