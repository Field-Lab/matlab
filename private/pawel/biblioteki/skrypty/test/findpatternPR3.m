%Sposob wywolywania NSReadClusterFile
%WaveformTypes=NS_ReadClusterFile(FileName,MovieNumber,PatternNumber,SD(1));

DataPath= 'D:\Home\Rydygier\Neuro\files';

ChannelNumber = 307;

PatternNumber = 30;
MovieNumber = 117;

MaxMovieNumber = 151;

FileName = [DataPath filesep 'ClusterFile_001_ID4546'];

%Define, which movie numbers to choose
FirstMovie = 91;
MovieStep = 5;
MovieNumbers = [145 147 149 151 153];


%Plotting data figures for movies defined in MovieNumbers
figure(1);
for i = 1:5
[DataTraces1,ArtifactDataTraces,DataChannels]=NS_ReadPreprocessedData(DataPath,DataPath,0,PatternNumber,MovieNumbers(i),0,0); 
DataTraces0=DataTraces1(1:100,:,:);
SD=size(DataTraces0);
ChannelTraces = DataTraces0(:,ChannelNumber,:);
ChannelTraces2D = reshape(ChannelTraces,SD(1),SD(3));
subplot(5,1,i), h= plot(ChannelTraces2D');

WaveformTypes=NS_ReadClusterFile(FileName,MovieNumbers(i),PatternNumber,SD(1));
artifactsIndex = find(WaveformTypes==1);
spikeIndex = find(WaveformTypes==2);
set(h(artifactsIndex),'Color','Black');
set(h(spikeIndex),'Color','Red');

h1=gca;
set(h1,'YLim',[-700 -200]);
grid on
end

%Plotting efficiency as a function of movie number
Spikes = zeros((MaxMovieNumber+1)/2,1);
k = 1;

for i=1:2:MaxMovieNumber
    WaveformTypes=NS_ReadClusterFile(FileName,i,PatternNumber,SD(1));
    Spikes(k,1) = i;
    Spikes(k,2) = sum(WaveformTypes==2)+sum(WaveformTypes==4);% || WaveformTypes==4);
    k = k+1;
end

AllMovies=[1:2:151];
Amps=1.05.^AllMovies;
figure(22);
plot(AllMovies,Spikes(:,2),'bd');
grid on;
xlabel('Movie NUmber');
ylabel('Stimulation Efficiency');

Effic = Spikes(:,2);

%Plotting stimulation efficiency as a function of movie number
PatternNumbers = [22 27 26 30 58 61];

figure(35);
for j=1:1
    %Data reading
    %[DataTraces1,ArtifactDataTraces,DataChannels]=NS_ReadPreprocessedData(DataPath,DataPath,0,PatternNumbers(i),MovieNumbers(i),0,0);
    k = 1;
    for i=1:2:MaxMovieNumber
        WaveformTypes=NS_ReadClusterFile(FileName,i,PatternNumbers(j),SD(1));
        Spikes(k,1) = i;
        Spikes(k,2) = sum(WaveformTypes==2)+sum(WaveformTypes==4);% || WaveformTypes==4);
        k = k+1;
    end
    Effic = Spikes(:,2);
    AllMovies = Spikes(:,1);
    subplot(6,1,j);
                   
    fo_ = fitoptions('method','NonLinearLeastSquares','Normalize','off','Robust','On');
    st_ = [0.71776047526892917 0.67238073502381523 ];
    set(fo_,'Startpoint',st_);
    ft_ = fittype('100/(1+exp(-x*a+b))',...
        'dependent',{'y'},'independent',{'x'},...
        'coefficients',{'a', 'b'});

    % Fit this model using new data
    cf_ = fit(Spikes(:,1),Effic(:,1),ft_,fo_);
    
    cf_.a
    cf_.b

    %Or use coefficients from the original fit:
    if 0
        cv_ = { 13.069381966060135, 13.33137736595198};
        cf_ = cfit(ft_,cv_{:});
    end
    
    cf_.a
    cf_.b

    
    a(j) = cf_.a;
    b(j) = cf_.b;
    for l = 1:length(Spikes)
        fit_data(l,1) = 100/(1+exp(-Spikes(l,1)*a(j)+b(j)));
    end
    plot(Spikes(:,1),fit_data);
    axis([0 160 0 100]);

end

% for i=1:64
%     [StimChannels,Amplitudes]=NS_StimulatedChannels(DataPath,i,1,[1:512],NS_GlobalConstants);
%     if find(StimChannels==299)
%         pattern=i
%     end
% end