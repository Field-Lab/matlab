function [OutputEvents,LatencyCriterion,HalfAmpCrossing]=NS512_FindSpikes(DataTraces,Channels,ThresholdNegative,ThresholdPositive,SpikeWidthRange,LatencyRange);
%- DataTraces - traces after artifact subtraction, by default - all the
%channels of the array (61 or 512)
%- Channels = on which channels spikes should be found
%ThresholdNegative,ThresholdPositive - obvious; the ThresholdNegative
%should be given WITH the "-"
%- SpikeWidthRange typically from 5 to 20: the width is between crossing of
%half of negative spike amplitude, and the first following crossing of zero
%- LatencyRange - the acceptable timing range for the crossing of half of
%negative spike amplitude (not threshold!!)

InputDataTracesSize = size(DataTraces); %Saving original input data size

DataTraces = DataTraces(:,Channels,:);
TracesSize = size(DataTraces); % no or ftraces, no of channels, no of samples 

EventsUnderNeg=DataTraces<=ThresholdNegative;
NegativeThresholdTraces=sign(sum(EventsUnderNeg,3)); %which traces have samples below negative threshold (2D array)

EventsOverPos = DataTraces>=ThresholdPositive;
PositiveThresholdTraces=sign(sum(EventsOverPos,3)); %which traces have samples above positive threshold (2D array)

BothThresholdTraces=NegativeThresholdTraces.*PositiveThresholdTraces;
%sum(BothThresholdTraces)

[TracesWithBothThresholds,ChannelsWithBothThresholds]=find(BothThresholdTraces==1);

SpikeWidthCriterion=zeros(TracesSize(1),TracesSize(2));
LatencyCriterion=zeros(TracesSize(1),TracesSize(2));
HalfAmpCrossing=LatencyCriterion;
for i=1:length(TracesWithBothThresholds)
    Tr=TracesWithBothThresholds(i);
    Ch=ChannelsWithBothThresholds(i);
    Trace=reshape(DataTraces(Tr,Ch,:),1,TracesSize(3));
    [a1,a2]=NS512_IsThisSpike(Trace,SpikeWidthRange); 
    SpikeWidthCriterion(Tr,Ch)=a1;
    HalfAmpCrossing(Tr,Ch)=a2;
    if a2>LatencyRange(1) & a2<LatencyRange(2)
        LatencyCriterion(Tr,Ch)=1;
    end        
end
OutputEvents=BothThresholdTraces.*SpikeWidthCriterion.*LatencyCriterion;