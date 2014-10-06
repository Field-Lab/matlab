% ten skrypt zaklada, ze istnieja juz opodwiednie pliki typu "ClusterFile",
% wygenerowane wczesniej przy (typowo recznej) analizie "podejrzanych"
% przypadkow, ktore jeszcze wczesniej wygenerowal skrypt
% "NS512_PatternAnalysis_2013_01_24" (moze byc z inna data)
%Procedura dla tego skryptu:
%1) Dla danego "podejrzanego" eventu wczytaj informacje z "CrosstalkFile"
%aby wiedziec o ktora stymulowana elektrode, oraz ktory neuron z
%potencjalnym crosstalkiem chodzi
%2) 
NS_GlobalConstants=NS_GenerateGlobalConstants(500);
Movies=[1:2:63];
EventNumber=82;
RawDataPath='F:\analiza\retina\2012-09-27-4\files\scan_new';
CrosstalkFilePath='C:\pawel\nauka\analiza\retina\2012-09-27-4\files\CrosstalkInfo';
ClusterFilePathCore='C:\pawel\nauka\analiza\retina\2012-09-27-4\files\ClusterFiles\ClusterFile_data000_ID'; %sciezka musi byc dopelniona numerem ID neuronu
EIfilePath='C:\pawel\nauka\analiza\retina\2012-09-27-4\files\EIs';
DiffEIsPath='C:\pawel\nauka\analiza\retina\2012-09-27-4\files\DiffEIs_try2';

CrosstalkInfo=NS512_ReadCrosstalkFile(CrosstalkFilePath,EventNumber)
Pattern=CrosstalkInfo(1)
Movie=CrosstalkInfo(2)
NeuronID=CrosstalkInfo(6)
PrimaryRecEl=CrosstalkInfo(7)

ClusterFilePath=[ClusterFilePathCore num2str(NeuronID)]
ClusterIndexes=NS_ReadClusterFileAll(ClusterFilePath)
eff=sum(ClusterIndexes-1,3);
EffVsAmp=eff(Movies,Pattern);
MoviesToAnalyze=Movies(find(EffVsAmp));
L=140;
t=[1:L];

electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(500); %define the electrode map - must be different than 1 for silicon probes
ChannelsPlot=electrodeMap.getAdjacentsTo(PrimaryRecEl,1)';

EI=NS512Read_EI_File(EIfilePath,NeuronID,ChannelsPlot);
SEI=size(EI);
for i=1:SEI(1)
    EI(i,:)=EI(i,:)-EI(i,1);
end
FigureProperties=struct('FigureNumber',120,'Subplot',[2 3 3],'TimeRange',[0 L],'AmplitudeRange',[min(EI(1,:))*1.5 max(EI(1,:))*1.5],'FontSize',16,'Colors',['k' 'r' 'b' 'm' 'g' 'c' 'y'],'LineWidth',2,'YLabel','Signal [mV]');

for i=1:length(MoviesToAnalyze)-1
    Movie1=MoviesToAnalyze(i);
    Movie2=MoviesToAnalyze(i+1);
    
    Movie1Indexes=find(ClusterIndexes(Movie1,Pattern,:)==1);
    Movie2Indexes=find(ClusterIndexes(Movie2,Pattern,:)==1);
    
    [DataTraces1,ArtifactDataTraces,Channels]=NS_ReadPreprocessedData(RawDataPath,RawDataPath,0,Pattern,Movie1,1500,Movie1Indexes);
    [DataTraces2,ArtifactDataTraces,Channels]=NS_ReadPreprocessedData(RawDataPath,RawDataPath,0,Pattern,Movie2,1500,Movie2Indexes);   
    
    SDT1=size(DataTraces1);
    p1=DataTraces1(:,ChannelsPlot,:);
    SDT2=size(DataTraces2);
    p2=DataTraces2(:,ChannelsPlot,:);
    DiffEI=mean(p2)-mean(p1);
    DiffEIAmp=max(DiffEI,[],3)-min(DiffEI,[],3);
    EIAmp=max(EI,[],2)-min(EI,[],2)    
    if (DiffEIAmp(1)>0.7*EIAmp(1))
        DataTraces=zeros(2,length(ChannelsPlot),L);
        DataTraces(1,:,:)=DiffEI;
        DataTraces(2,:,:)=EI;
        WaveformTypes=[1 0];
        y=NS_PlotClustersOfSignaturesOnArrayLayout(DataTraces,ChannelsPlot,WaveformTypes,500,FigureProperties,NS_GlobalConstants);
        
        %DiffEIFileName=[DiffEIsPath '\p' num2str(Pattern) '_m' num2str(Movie2) '_ID' num2str(NeuronID) '_rec' num2str(PrimaryRecEl) '_e' num2str(EventNumber)]                
        %h=gcf;
        %set(h,'PaperUnits','inches');
        %set(h,'PaperSize',[16 9]);
        %set(h,'PaperPosition',[0 0 16 9]); 
        %print(h, '-dtiff', '-r120', DiffEIFileName);        
    end
end