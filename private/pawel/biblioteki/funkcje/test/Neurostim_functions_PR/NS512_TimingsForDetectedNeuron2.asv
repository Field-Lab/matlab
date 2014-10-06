function [CorrectedTraces,EI,TimingsForSpikesMinimumsCorrected]=NS512_TimingsForDetectedNeuron2(TracesWithSpikes,ChannelsPlot);
%This function calculates 

% 1. Find spikes timings - index of most negative sample for each spike
% (simplified version now!! PH, 2010-06-07)

ChannelTraces=TracesWithSpikes(:,ChannelsPlot(1),:);
SCT=size(ChannelTraces);
ChannelTraces2D=reshape(ChannelTraces,SCT(1),SCT(3));
TimingsForSpikesMinimums=FindUnifiedSpikes_PR(ChannelTraces2D,-15); %Wektor wskaznikow (numery probek dla ktorych nastapilo przekroczenie progu) - zmienic na: przekroczona zostala polowa amplitudy

% 3. Find spikes with timing within specified range, and correct timing if
% some spikes are too early
[SpikesWithinTimeBrackets,Indexes]=NS512_SpikesWithinTimingBrackets(TimingsForSpikesMinimums,1.5,1.5); % 1 - given spike is within the range, 0 - it is not
%
TracesWithSpikesWithinBrackets=TracesWithSpikes(Indexes); %numbers of traces are with spikes that are within brackets :)
TimingsForSpikesMinimumsCorrected=NS512_CorrectSpikesTimings(TimingsForSpikesMinimums(Indexes),Indexes); %correct the spikes timings if they were too soon

% 4. Find traces after jitter correction and the EI
DataForEI=TracesWithSpikes(Indexes,ChannelsPlot,:);
[CorrectedTraces,EI,Noise]=NS512_CorrectSignaturesForTiming(DataForEI,TimingsForSpikesMinimumsCorrected);