function [timeLog,psth] = psth_calc(trialData,binSize,windowing)
numTrials=size(trialData,1);
trialLen=size(trialData,2);
timeLog=[];

if(strcmp(windowing,'overlap'))
psthLen=trialLen-binSize+1;
    psth=zeros(1,psthLen);
for itime=1:psthLen
    psth(itime)=sum(sum(trialData(:,itime:itime+binSize-1)));
    timeLog=[timeLog;itime+(binSize/2)];
end

end

if(strcmp(windowing,'nonoverlap'))
psthLen = floor(trialLen / binSize);
icnt=0;
psth=zeros(1,psthLen);

for itime=1:binSize:trialLen-binSize+1
icnt=icnt+1;
    psth(icnt)= sum(sum(trialData(:,itime:itime+binSize-1)));
timeLog=[timeLog;itime+(binSize/2)];
end

end