NS_GlobalConstants=NS_GenerateGlobalConstants(512);

FilePath='E:\pawel\analysis\cultures\2009-11-20-0\2009-11-20-0\data003e\';
FilePath='D:\analysis\cultures\2009-11-20-0\data003e\';

%Amplitudes:
for m=1:0%55
    [StimChannels,Amplitudes]=NS_StimulatedChannels(FilePath,1,m,[1:512],NS_GlobalConstants);
    A(m)=max(max(abs(Amplitudes)));
end

NeuronNumber=5;

figure(1);
for i=1:30
    FileName=['ClusterFile_003_n' num2str(i)];
    FullPath=[FilePath FileName];
    ClusterIndexes=NS_ReadClusterFileAll(FullPath);
    CI=ClusterIndexes(:,:,1:100);
    h=NS512_PlotEffVsAmp(CI);
    subplot(5,6,i);
end