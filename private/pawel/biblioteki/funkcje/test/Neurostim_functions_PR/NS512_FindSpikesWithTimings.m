function [TracesWithoutArtifact,Artifact,Events,HalfAmpCrossingIndexes] = NS512_FindSpikesWithTimings(Traces,Channels,ArtifactThresholdNumber,SpikesNumberThreshold);

% 1. Artifact estimation
[Artifact]=ArtifactEstimation_PR(Traces,ArtifactThresholdNumber); %size of Artifact is 1 x 512 x 70

% 2. Artifact subtraction 
TracesWithoutArtifact=NS512_SubtractArtifact(Traces,Artifact); %size of TracesWithoutArtifact is 100 x 512 x 70

% 3. Spike detection
[Events,LatencyCriterion,HalfAmpCrossingIndexes]=NS512_FindSpikes(TracesWithoutArtifact,Channels,-25,3,[5 20],[3 40]);