function [TracesWithoutArtifact,Artifact,Events,ChannelsWithSpikes,SpikesTimings] = FindTracesClasses_PR2(Traces,Channels,ArtifactThresholdNumber,SpikesNumberThreshold);

% 1. Artifact estimation
[Artifact]=ArtifactEstimation_PR(Traces,ArtifactThresholdNumber); %size of Artifact is 1 x 512 x 70

% 2. Artifact subtraction 
TracesWithoutArtifact=NS512_SubtractArtifact(Traces,Artifact); %size of TracesWithoutArtifact is 100 x 512 x 70

% 3. Spike detection
%[Events,~] = FindSpikes2_PR(TracesWithoutArtifact,Channels, -25, 3, 3, 20); %add possibility to specify list of channels for spike finding (PH, 2010-06-01)
[Events,LatencyCriterion,HalfAmpCrossingIndexes]=NS512_FindSpikes(TracesWithoutArtifact,Channels,-25,3,[5 20],[3 40]);

% 4. Find channels with more than 50 spikes
EventsPerChannel=sum(Events); %How many spikes on each electrode, size - 1 x 512
ChannelsWithSpikes0=find(EventsPerChannel > SpikesNumberThreshold); % ktore elektrody zarejestrowaly wiecej niz 50 spikow
ChannelsWithSpikes=Channels(ChannelsWithSpikes0);


% 5. Find timings for spikes on specified channels
%SpikesTimings=NS512_FindSpikesTimings(TracesWithoutArtifact,Events,ChannelsWithSpikes); %separately for each neuron = each element of the ChannelsWithSpikes array
%asdasd=size(SpikesTimings)