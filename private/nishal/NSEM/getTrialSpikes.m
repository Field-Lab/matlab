function [spikeTrial,spikeTrialBinned]=getTrialSpikes(dataTrial,spikeThreshold,thrMethod,binSize)

if(thrMethod=='max')
spikeTrial=double((dataTrial>spikeThreshold));
else if(thrMethod=='min')
spikeTrial=double((dataTrial<spikeThreshold));
    end
end
%Remove duplicates , etc ?? 

spikeTrialBinned=zeros(size(dataTrial,1),size(dataTrial,2)/binSize);

for iTrial=1:size(dataTrial,1)
    icnt=1;
    for iSample=1:binSize:size(dataTrial,2)-binSize     
        icnt=icnt+1;
        spikeTrialBinned(iTrial,icnt)= sum(spikeTrial(iTrial,iSample:iSample+binSize))>0;
    end
end

end