clear;
NS_GlobalConstants=NS_GenerateGlobalConstants(61);
p='E:\2008-08-26-0\';
FileName='data005\data005000.bin';
paramsFile=edu.ucsc.neurobiology.vision.io.ParametersFile('E:\analysis\2008-08-26-0\data005\data005.params');
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile('E:\analysis\2008-08-26-0\data005\data005.neurons');

%p='D:\2008-08-27-4\';
%FileName='data004\data004000.bin';
%paramsFile=edu.ucsc.neurobiology.vision.io.ParametersFile('D:\2008-08-27-4\analysis\data004\data004.params');
%neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile('D:\2008-08-27-4\analysis\data004\data004.neurons');

idList = neuronFile.getIDList();

NeuronID=602;
spikeTimes = neuronFile.getSpikeTimes(NeuronID)';
Timings1=spikeTimes(1,1:min(1000,length(spikeTimes)))-7;
Channels=[1:64];
CenterChannel=41;
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(1);
Radius=1;
ChannelsPlot=electrodeMap.getAdjacentsTo(CenterChannel,Radius)';
%Channels=[44 2];
[RAWtraces,signal]=NS_AverageTraces([p FileName],Timings1-1,Channels,[-3 57],NS_GlobalConstants);
signal=signal';
ss=size(signal);
for i=1:ss(1)
    signal(i,:)=signal(i,:)-mean([signal(i,ss(2)) signal(i,1)]);
end
%signal=signal*2;
s(1,:,:)=signal(ChannelsPlot,:);

FigureProperties=struct('FigureNumber',104,'Subplot',[2 3 3],'TimeRange',[1 55],'AmplitudeRange',[-300 200],'FontSize',24,'Colors',['k' 'r' 'b' 'k' 'g' 'm' 'c'],'LineWidth',2,'YLabel','input signal [\muV]');
y=NS_PlotClustersOfSignaturesOnArrayLayout(s/0.44,ChannelsPlot,[1],1,FigureProperties,NS_GlobalConstants);

break;

TracesNew=NS_NormalizeSignatures(s,100);
figure(106);
clf
el=[49 62];
for i=1:0 %length(el)
    figure(106);
    plot(signal(el(i),:));
    hold on;
end
grid on

ss=size(s);
sig=reshape(s,ss(2),ss(3));
FigureProperties.FigureNumber=105;
Channels=[1:64];
M=NS_SaveMovieFromSignature(signal,Channels,CenterChannel,1,FigureProperties,NS_GlobalConstants);
cd C:\praca\analiza\report_2008_10_16\movies;
movie2avi(M,['Neuron' num2str(NeuronID)],'fps',15,'quality',100);