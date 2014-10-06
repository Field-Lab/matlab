function MainClusterIndexes=NS512_FindMainClusterIndexesForSuspectedCrosstalkEvents(Data,MatchingThreshold);
%Zakladamy ze Data to dane z JEDNEGO KANALU (primary recording electrode
%dla danego eventu) - 2013/03/20
%output: indexes for all the traces that belong to the main cluster. For
%now (2013/03/20) it is just one trace!
SDT=size(Data);

UnMatched=zeros(SDT(1),1);
for i=1:SDT(1)    
    for j=1:SDT(1)
        difference=Data(j,:)-Data(i,:);
        %DataPrim(j,:)=Data(j,:)-Data(i,:);
        ThresholdCross=find((abs(difference)>MatchingThreshold));
        if ThresholdCross
            UnMatched(i,1)=UnMatched(i,1)+1;
        end                            
    end   
end
[value,index]=min(UnMatched)
UnMatched
MainClusterIndexes=[];
for j=1:SDT(1)
    difference=Data(j,:)-Data(index,:);    
    ThresholdCross=find((abs(difference)>MatchingThreshold));
    if length(ThresholdCross)==0
        MainClusterIndexes=[MainClusterIndexes j];
        %UnMatched(i,1)=UnMatched(i,1)+1;
    end                            
end   


%difference=Data(j,:)-Data(i,:);
%MainClusterIndexes=index;