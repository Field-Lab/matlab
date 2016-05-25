electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(500);

DataPath='D:\Home\Pawel\analysis\slices\2013\2013-12-12-3-PH\data001preproc'
ArtifactDataPath=DataPath;
PatternsToAnalyze=[316];
MovieStart=13;

figure(11);
clf;
Pattern=211;
PatternX=electrodeMap.getXPosition(Pattern);
PatternY=electrodeMap.getYPosition(Pattern);
for m=MovieStart:18:(MovieStart+18*24)
    clf
    [DataTraces,ArtifactDataTraces,Channels]=NS_ReadPreprocessedData(DataPath,ArtifactDataPath,0,Pattern,m,1500,0);
    Events=NS512_DetectSpikes(DataTraces,25,10);
    [a,b]=find(Events>0);
    n=hist(b,[1:512]);
    ElectrodesWithSpikes=find(n>9)
    clear X;
    clear Y;
    for i=1:length(ElectrodesWithSpikes)
        X(i)=electrodeMap.getXPosition(ElectrodesWithSpikes(i));
        Y(i)=electrodeMap.getYPosition(ElectrodesWithSpikes(i));        
    end
    subplot('position',[0.5 0.5 0.05 0.05]);    
    if ElectrodesWithSpikes
        h1=plot(X,Y,'ro');
        set(h1,'MarkerFaceColor','r')
        set(h1,'MarkerSize',4)
    end
    hold on
    h1=plot(PatternX,PatternY,'gd');    
    set(h1,'MarkerFaceColor','g')
    axis([-1000 1000 -500 500])
    h=gca
    set(h,'XTickLabel','')
    set(h,'YTickLabel','')
    %pause(0.1)
    %FullName=['D:\Home\Pawel\analysis\slices\2013\2013-12-12-3-PH\data001_Matlab\SpikeAnalysis\proba1.tif'];   
    FullName=['D:\Home\Pawel\analysis\slices\2013\2013-12-12-3-PH\data001_Matlab\SpikeAnalysis\proba' num2str(m) '.tif']; 
    h=gcf;
    set(h,'PaperUnits','inches');
    set(h,'PaperSize',[16 9]);
    set(h,'PaperPosition',[0 0 16 9]); 
    print(h, '-dtiff', '-r120', FullName);
end
FullName=['D:\Home\Pawel\analysis\slices\2013\2013-12-12-3-PH\data001_Matlab\SpikeAnalysis\proba' num2str(m) '.tif'];            
h=gcf;
set(h,'PaperUnits','inches');
set(h,'PaperSize',[16 9]);
set(h,'PaperPosition',[0 0 16 9]); 
%set(h,'PaperPosition',[0 0 10.7 13.5]); 
print(h, '-dtiff', '-r120', FullName);