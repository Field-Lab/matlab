clear

NS_GlobalConstants=NS_GenerateGlobalConstants(500);

neuronFile=edu.ucsc.neurobiology.vision.io.NeuronFile('D:\Home\Pawel\analysis\slices\2013\2013-12-12-3-PH\data005\data005.neurons');
rawFile=edu.ucsc.neurobiology.vision.io.RawDataFile('I:\analysis\slices\2013-12-12-3-PH\data005');


idList=[96 216 217 289];
Colors={'b' 'r' 'g' 'k'};
figure(10)
clf

CenterChannel=15;
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(500); %define the electrode map - must be different than 1 for silicon probes
Radius=1;
ChannelsPlot=electrodeMap.getAdjacentsTo(CenterChannel,Radius)';

N=100;
L=100;
spikes=zeros(N*length(idList),7,L);
WaveformTypes=zeros(1,length(idList)*N);
for n=1:length(idList)
    NeuronID=idList(n)
    spikeTimes = neuronFile.getSpikeTimes(NeuronID)';        
    
    for i=1:N
        t=spikeTimes(i+5);
        d0=rawFile.getData(t-40,L)';
        d1=d0(ChannelsPlot+1,:);
        spikes(i+(n-1)*N,:,:)=d1-(n-1)*120;
    end
    WaveformTypes(1,(n-1)*N+1:n*N)=n;
    %h=plot(spikes'-(n-1)*100);
    %set(h,'Color',Colors{n});
    %hold on
end

FigureProperties=struct('FigureNumber',12,'Subplot',[2 3 3],'TimeRange',[1 100],'AmplitudeRange',[-500 50],'FontSize',14,'Colors',['k' 'r' 'b' 'k' 'g' 'm' 'c'],'LineWidth',1,'YLabel','input signal [\muV]');
y=NS_PlotClustersOfSignaturesOnArrayLayout(spikes,ChannelsPlot,WaveformTypes,500,FigureProperties,NS_GlobalConstants);

FullImageName=['C:\home\Pawel\nauka\analiza\report_2014_01_23\RawTraces2.tif'];
h=gcf;
set(h,'PaperUnits','inches');
set(h,'PaperSize',[10 10]);
set(h,'PaperPosition',[0 0 10 10]); 
print(h, '-dtiff', '-r240', FullImageName);