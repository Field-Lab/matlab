function [Events,Artifact]=NS512FindResponsesToLocalStimulation(DataTraces,ThresholdNeg,ThresholdPos,N);
%- Events - 2-D array which says which traces include spikes, independently
%for each channel
%Traces - NumberOfTraces x NumberOfChannels x NumberOfSamples
%N - how many traces are used for artifact estimation
[Traces,Artifact,c]=NS512_SubtractLocalArtifact(DataTraces,N);
Events=NS512_DetectSpikes(Traces,ThresholdNeg,ThresholdPos);