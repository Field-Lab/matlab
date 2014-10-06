function Artifact=NS512_ArtifactTracesInSomaticResponses(Traces,Indexes);
%Traces - 2-D array: repetitions,samples
%Indexes - result of sorting: first traces are spikes, last traces are
%artifacts only, but how many are spikes and how many artifacts? This
%function finds how many last traces are artifacts.

ST=size(Traces);
for i=1:ST(1)
    