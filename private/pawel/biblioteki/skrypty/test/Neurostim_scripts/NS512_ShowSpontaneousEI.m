%D:\analysis\cultures\2009-11-20-0\data002
NS_GlobalConstants=NS_GenerateGlobalConstants(512);
FileName='H:\2009-11-27-0\data000';
%paramsFile=edu.ucsc.neurobiology.vision.io.ParametersFile('C:\home\pawel\praca\analysis\2009-11-20-0\data002\data002.params');
paramsFile=edu.ucsc.neurobiology.vision.io.ParametersFile('F:\analysis\retina\2009-11-27-0\data000\data000000\data000000.params');
%neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile('C:\home\pawel\praca\analysis\2009-11-20-0\data002\data002.neurons');
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile('F:\analysis\retina\2009-11-27-0\data000\data000000\data000000.neurons');
%eiFile=edu.ucsc.neurobiology.vision.io.EIFile('D:\analysis\cultures\2009-11-20-0\data002\data002.ei');
idList = neuronFile.getIDList();
NeuronID=3631;
spikeTimes = neuronFile.getSpikeTimes(NeuronID)';
N=500;
L=100;
Timings1=spikeTimes(1,1:min(N,length(spikeTimes)))-20;
%break;
Data=zeros(513,L);
rawFile=edu.ucsc.neurobiology.vision.io.RawDataFile(FileName); 
for i=1:length(Timings1)
    RawData=double(rawFile.getData(Timings1(i),L)');    
    Data=Data+RawData;
end
m=mean(Data')';
for i=1:513
    Data(i,:)=Data(i,:)-m(i);
end
%clear T;

T(1,:,:)=Data(2:513,:)/N;
FigureProperties=struct('FigureNumber',3,'TimeRange',[0 80],'AmplitudeRange',[-60 40],'FontSize',16,'Colors',['k' 'r' 'b' 'y'],'XTick',[0 0.5 1 1.5 2 2.5 3 3.5],'YTick',[-60:20:40],'LineWidth',2,'YLabel','signal');
%y=NS512_PlotClustersOfSignaturesOnArrayLayoutWithMarks(T,[425:512],[1],500
%,FigureProperties,NS_GlobalConstants,[484 493 459 434 504 499]);