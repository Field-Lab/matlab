function [TracesWithoutArtifact,Artifact,Events,ChannelsWithSpikes, TimeCorelSpikes] = FindTracesClasses_PR(Traces,Channels,ArtifactThresholdNumber,SpikesNumberThreshold );

% 1. Artifact estimation
[ Artifact ] = ArtifactEstimation_PR(Traces,ArtifactThresholdNumber); %size of Artifact is 1 x 512 x 70
figure(103)
size(Artifact)
plot(reshape(Artifact(1,130,:),1,70));

% 2. Artifact subtraction 
TracesWithoutArtifact = NS512_SubtractArtifact(Traces,Artifact); %size of TracesWithoutArtifact is 100 x 512 x 70

% 3. Spike detection
[Events,~] = FindSpikes2_PR(TracesWithoutArtifact,Channels, -25, 3, 3, 20); %add possibility to specify list of channels for spike finding (PH, 2010-06-01)

% 4. Find channels with more than 50 spikes
EventsPerChannel=sum(Events); %How many spikes on each electrode, size - 1 x 512
ChannelsWithSpikes=find(EventsPerChannel > SpikesNumberThreshold); % ktore elektrody zarejestrowaly wiecej niz 50 spikow

% 5. Find timings for spikes for specified channels
SpikesTimings=NS512_FindSpikesTimings(TracesWithoutArtifact,Events,ChannelsWithSpikes);

% 6. Find spikes with timing within specified range
TimeCorelSpikes = FindTimeCorellatedSpikes_PR(SpikesTimings, 1.5, 1.5); % 1 - given spike is within the range, 0 - it is not

% 7. Find indexes for electrophysiological imaging
Indexes=NS512_IndexesForSignatures(SpikesTimings,TimeCorelSpikes);