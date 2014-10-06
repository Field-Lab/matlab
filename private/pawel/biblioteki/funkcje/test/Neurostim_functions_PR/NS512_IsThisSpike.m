function [Answer,HalfAmpCrossing]=NS512_IsThisSpike(Trace,SpikeWidthRange);
%For a waveform that can be a spike or not, the function checks whether the
%delay between first negative crossing of the half of the negative
%amplitude, and the first positive crossing of zero, is within defined
%range.
%Answer = 0 or 1
%HalfAmpCrossing = when the first crossing of the half of the negative amplitude
%happens

SamplesUnderHalfAmp=Trace<min(Trace)/2;
HalfAmpCrossing=find(diff(SamplesUnderHalfAmp)==1,1)+1; %can happen more than once! so take the first index
if length(HalfAmpCrossing)==0
    Answer=0;
    HalfAmpCrossing=0;
    return;
end

SamplesAboveZero=Trace>0;
SamplesAboveZero=Trace(HalfAmpCrossing:length(Trace))>0;

ZeroCrossing=find(diff(SamplesAboveZero)==1,1)+1;

SpikeWidth=ZeroCrossing;
if SpikeWidth>SpikeWidthRange(1) & SpikeWidth<SpikeWidthRange(2)
    Answer=1;
else
    Answer=0;
end