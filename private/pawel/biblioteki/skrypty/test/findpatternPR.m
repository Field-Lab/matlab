%Sposob wywolywania NSReadClusterFile
%WaveformTypes=NS_ReadClusterFile(FileName,MovieNumber,PatternNumber,SD(1));

DataPath= 'D:\Home\Rydygier\Neuro\files';

ChannelNumber = 307;

PatternNumber = 26;
MovieNumber = 117;

MaxMovieNumber = 151;

FileName = [DataPath filesep 'ClusterFile_001_ID4546'];


%Define, which movie numbers to choose
FirstMovie = 91;
MovieStep = 5;
MovieNumbers = [147 149 151 153];
%Plotting data figures for movies defined in MovieNumbers
figure(1);
for i = 1:4
[DataTraces1,ArtifactDataTraces,DataChannels]=NS_ReadPreprocessedData(DataPath,DataPath,0,PatternNumber,MovieNumbers(i),0,0); 
DataTraces0=DataTraces1(1:100,:,:);
SD=size(DataTraces0);
ChannelTraces = DataTraces0(:,ChannelNumber,:);
ChannelTraces2D = reshape(ChannelTraces,SD(1),SD(3));
subplot(2,2,i), h= plot(ChannelTraces2D');

WaveformTypes=NS_ReadClusterFile(FileName,MovieNumbers(i),PatternNumber,SD(1));
artifactsIndex = find(WaveformTypes==1);
spikeIndex = find(WaveformTypes==2);
set(h(artifactsIndex),'Color','Black');
set(h(spikeIndex),'Color','Red');
           
end

%Plotting efficiency as a function of movie number
Spikes = zeros((MaxMovieNumber+1)/2,1);
k = 1;
for i=1:2:MaxMovieNumber
    WaveformTypes=NS_ReadClusterFile(FileName,i,PatternNumber,SD(1));
    Spikes(k,1) = i;
    Spikes(k,2) = sum(WaveformTypes==2);
    k = k+1;
end

figure(22);
plot(Spikes(:,1),Spikes(:,2),'bd');
grid on;
xlabel('Movie NUmber');
ylabel('Stimulation Efficiency');

% for i=1:64
%     [StimChannels,Amplitudes]=NS_StimulatedChannels(DataPath,i,1,[1:512],NS_GlobalConstants);
%     if find(StimChannels==299)
%         pattern=i
%     end
% end