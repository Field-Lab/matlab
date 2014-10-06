Pattern=282;
Movie=53;
Electrode=105;

StimulationDataPath='F:\analiza\retina\2012-09-27-4\files\scan_new';
[DataTraces0,ArtifactDataTraces,DataChannels]=NS_ReadPreprocessedData(StimulationDataPath,StimulationDataPath,0,Pattern,Movie,0,0);   %read data                 
Data=reshape(DataTraces0(1:50,Electrode,:),50,140);

prog=10;
SDT=size(Data);

figure(2);
UnMatched=zeros(SDT(1),1);
for i=1:SDT(1)    
    for j=1:SDT(1)
        difference=Data(j,:)-Data(i,:);
        %DataPrim(j,:)=Data(j,:)-Data(i,:);
        ThresholdCross=find((abs(difference)>prog));
        if ThresholdCross
            UnMatched(i,1)=UnMatched(i,1)+1;
        end                            
    end
    [value,index]=min(UnMatched)
    BestTraceIndex=index;
end
figure(1)
clf
plot(Data');
hold on
h=plot(Data(index,:));
set(h,'LineWidth',5)
h=gca;
set(h,'XLim',[0 100])
set(h,'YLim',[-440 -260])