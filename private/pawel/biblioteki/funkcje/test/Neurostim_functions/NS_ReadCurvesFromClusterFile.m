function [Amplitudes,Curves]=NS_ReadCurvesFromClusterFile(DataPath,FileName,CenterChannel,Patterns,Movies,NS_GlobalConstants);

electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(1);
Channels=electrodeMap.getAdjacentsTo(CenterChannel,1)';
for i=1:length(Movies)
    FullName_status=[DataPath filesep 'status_m' num2str(Movies(i))];
    load(FullName_status);
    FullName_pattern=[DataPath filesep 'pattern' num2str(Patterns(1)) '_m' num2str(Movies(i))];
    load(FullName_pattern);
    [patternTimes,traces] = plotStimulusTraces(Pattern, Status, Channels, [0 50], NS_GlobalConstants);
    Amplitudes(i)=max(max(max(abs(traces))));
end
WaveformTypesAll=NS_ReadClusterFileAll(FileName);
Curves=zeros(6,length(Movies));

for i=1:length(Patterns)
    for j=1:length(Movies)
        Patterns(i);
        Movies(j);
        WaveformTypes=reshape(WaveformTypesAll(Movies(j),Patterns(i),:),100,1);
        a=find(WaveformTypes==1); % it is assumed that in the cluster files, for each clustering the cluster number 1 cirrespond to "artifact only"
        Curves(i,j)=100-length(a);
    end
end