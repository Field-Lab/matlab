%Sposob wywolywania NSReadClusterFile
%WaveformTypes=NS_ReadClusterFile(FileName,MovieNumber,PatternNumber,SD(1));
NS_GlobalConstants=NS_GenerateGlobalConstants(512);
DataPath= 'D:\Home\Rydygier\Neuro\files';

ChannelNumber = 307;

PatternNumber = 26;
MovieNumber = 117;

MaxMovieNumber = 151;

FileName = [DataPath filesep 'ClusterFile_001_ID4546'];

%Define, which movie numbers to choose
FirstMovie = 91;
MovieStep = 5;
MovieNumbers = [145 147 149 151 153];
%MovieNumbers = [71 75 83 89 89];

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
set(h1,'YLim',[-700 200]);
grid on
end
%break;
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
for j=1:6
   k = 1;
    for i=1:2:MaxMovieNumber
        WaveformTypes=NS_ReadClusterFile(FileName,i,PatternNumbers(j),SD(1));
        Spikes(k,1) = i;
        Spikes(k,2) = sum(WaveformTypes==2)+sum(WaveformTypes==4);% || WaveformTypes==4);
        k = k+1;
    end
    Effic = Spikes(:,2);
    AllMovies = Spikes(:,1);
    
    Amplitudes =FindAmplitudes_PR( DataPath,1,[1:2:151],[1:512],NS_GlobalConstants );
    cf_=SigmoidalFit_PR(Amplitudes,Effic);
    figure(36);

%Dla ustawienia w poziomie
%    subplot(1,6,j);
%Dla ustawienia w pionie
     subplot(6,1,j);

    h=plot(Amplitudes, Spikes(:,2),'bd','Markersize',6);
    set(h,'LineWidth',0.1);
    hold on;
    h=plot(cf_);
    set(gca,'Xscale','log');
    set(h,'LineWidth',2);
    legend off;
    
    axis([0.1 4 0 100]);
    set(gca,'XTick',[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 2 3 4]);
    %set(gca,'XTickLabel',{'0.1','','0.3','','','','','','','1','','3',''});
    
    set(gca,'YTick',[0 25 50 75 100]);
    set(gca,'YTickLabel',{'0','','50','','100'});
    set(gca,'GridLineStyle','-');
    grid off;
%Dla ustawienia w pionie

    
     if j== 6
         xlabel(gca,'Amp. [\muA]','FontSize',12);
         ylabel(gca,'Eff. [%]','FontSize',12);
         set(gca,'XTickLabel',{'0.1','','0.3','','','','','','','1','','3',''});
     else
         xlabel(gca,'');
         ylabel(gca,'');
         set(gca,'XTickLabel',{'','','','','','','','','','','','',''});
     end
     set(gca,'FontSize',12);

%Dla ustawienia w poziomie
%{    
    if j==1
        ylabel(gca,'Eff.[%]','FontSize',16);
        xlabel(gca,'Amp.[\muA]','FontSize',16);
    else
        ylabel(gca,'');
        xlabel(gca,'');
    end
    set(gca,'FontSize',14);
%}    
    %figure(35);
    %a(j) = cf_.a;
    %A=a(j);
    %b(j) = cf_.b;
    %B=b(j);
    %for l = 1:length(Spikes)
    %    fit_data(l,1) = 100/(1+exp(-Spikes(l,1)*a(j)+b(j)));
    %    fit_data(l,1) = 100/(1+exp(-AllMovies*A+B));
    %end
    %fit_data=100./(1+exp(-AllMovies*A+B));
    %plot(Spikes(:,1),fit_data);
    %plot(fit_data);
    %axis([0 2 0 100]);
    hold off;
end

FullName = 'D:\Home\Rydygier\NEURO\analysis\retina\2009-11-27-0\data001\April2010\figures\toPaper\sigma_fit_vertical_log3';

hc = figure(36);
set(hc,'PaperUnits','inches');
        set(hc,'PaperSize',[2 9]);
        set(hc,'PaperPosition',[0 0 2 9]);
        
        print(hc, '-dtiff', '-r120', FullName);

% for i=1:64
%     [StimChannels,Amplitudes]=NS_StimulatedChannels(DataPath,i,1,[1:512],NS_GlobalConstants);
%     if find(StimChannels==299)
%         pattern=i
%     end
% end