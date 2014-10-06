StimulationDataPath='F:\analiza\retina\2012-09-27-4\files\scan_new';
CrosstalkFilePath='C:\pawel\nauka\analiza\retina\2012-09-27-4\files\CrosstalkInfo';
EventMax=50;

MatchingThreshold=20;

for i=[113]% 115 129 132]%[2:99 101:214]
    CrosstalkInfo=NS512_ReadCrosstalkFile(CrosstalkFilePath,i);
    Pattern=CrosstalkInfo(1);
    Movie=CrosstalkInfo(2);
    NeuronID=CrosstalkInfo(6);
    PrimaryRecEl=CrosstalkInfo(7);
    
    MovieMin=max(Movie-6,1);
    MovieMax=min(Movie+6,63);
    for M=MovieMin:2:MovieMax
        [DataTraces,ArtifactDataTraces,DataChannels]=NS_ReadPreprocessedData(StimulationDataPath,StimulationDataPath,0,Pattern,M,0,0);
        SDT=size(DataTraces);
        NumberOfEvents=min(SDT(1),EventMax);
        Data=reshape(DataTraces(1:NumberOfEvents,PrimaryRecEl,:),NumberOfEvents,SDT(3));
        MainClusterIndexes=NS512_FindMainClusterIndexesForSuspectedCrosstalkEvents(Data,MatchingThreshold);
        figure(2);
        clf
        plot(Data');
        hold on
        h=plot(Data(MainClusterIndexes,:)');
        set(h,'LineWidth',5)
        h=gca;
        set(h,'XLim',[0 100])
        set(h,'YLim',[-600 -200])
        
        h=gcf;
        FullName=['C:\pawel\nauka\analiza\retina\2012-09-27-4\files\AutomatedMainClusterFinding_try3\e' num2str(i) '_p' num2str(Pattern) '_m' num2str(M) '_rec' num2str(PrimaryRecEl)];
        set(h,'PaperUnits','inches');
        set(h,'PaperSize',[16 9]);
        set(h,'PaperPosition',[0 0 16 9]);  
        print(h, '-dtiff', '-r120', FullName);        
    end
end