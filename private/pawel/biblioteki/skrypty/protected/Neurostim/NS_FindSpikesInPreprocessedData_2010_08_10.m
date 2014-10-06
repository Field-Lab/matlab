ChipAddresses=[24:31];
NumberOfChannelsPerChip=64;
CurrentRanges=[0.066 0.266 1.07 4.25 16.9 67.1 264 1040];
Fs=20000;
%NS_GlobalConstants=struct('SamplingFrequency',Fs,'ChipAddresses',ChipAddresses,'NumberOfChannelsPerChip',NumberOfChannelsPerChip,'CurrentRanges',CurrentRanges);
ArrayID=500;

NS_GlobalConstants=NS_GenerateGlobalConstants(512);

WritePath='C:\home\pawel\2010\analysis\07_2010_Cultures\SpikesInStimData';
PrepDataPath='I:\analysis\2010-07-29-0\data';
ArtifactDataPath='I:\analysis\2010-07-29-0\artifacts';
ClusterFileName='I:\analysis\2010-07-29-0\artifacts';
GoodChannels=[1:512];
Movies=[1:8:209];

Patterns=NS512_OptimalElectrodeSequence();
PatternsForMovie=Patterns(1:64)
ArtifactSubtraction=1;
Responses=zeros(length(PatternsForMovie),512);
Noise=Responses;

tic
for PatternNumber=PatternsForMovie(1,1:4)
    PatternNumber
    [a,b]=NS512_ApplyLinearArtifactModel3(PrepDataPath,ArtifactDataPath,ArtifactSubtraction,PatternNumber,Movies,GoodChannels,0,0,ClusterFileName,NS_GlobalConstants);
    Responses(PatternNumber,:)=a;
    Noise(PatternNumber,:)=b;
end
toc
save Responses;
Save Noise;