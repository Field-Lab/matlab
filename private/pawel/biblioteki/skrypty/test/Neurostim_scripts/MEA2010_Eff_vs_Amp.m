NS_GlobalConstants=NS_GenerateGlobalConstants(512);
DataPath= 'D:\Home\Rydygier\Neuro\files';
ArrayID=500;

ChannelNumber = 242;
PatternNumber = 42;
MovieNumber = 117;
MaxMovieNumber = 152;

FileName = [DataPath filesep 'ClusterFile_001_el242_100us'];

%Define, which movie numbers to choose
FirstMovie = 91;
MovieStep = 5;

MovieNumbers = [32 38 44 50 58 66];

time=[-0.35:0.05:3.6]-0.1;

x=[0 1 1 3 3 5 5 7 7 8];
x=x*0.05-0.35;
y=[0 0 2 2 -3 -3 1 1 0 0];
Skalowanie=100;

electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(500); %define the electrode map - must be different than 1 for silicon probes
Radius=1;
ChannelsPlot=electrodeMap.getAdjacentsTo(ChannelNumber,Radius)';

%Plotting data figures for movies defined in MovieNumbers
figure(1);
%break;
for i=1:0 %length(MovieNumbers)
    Amplitudes =FindAmplitudes_PR( DataPath,1,MovieNumbers(i),[1:512],NS_GlobalConstants );
    CurrentAmplitude = Amplitudes(max(size(Amplitudes)));
    [DataTraces1,ArtifactDataTraces,DataChannels]=NS_ReadPreprocessedData(DataPath,DataPath,0,PatternNumber,MovieNumbers(i),0,0); 
    DataTraces0=DataTraces1(1:100,:,:)+350;
    SD=size(DataTraces0);
    Traces = DataTraces0(:,ChannelsPlot,:);    
    
    WaveformTypes=NS_ReadClusterFile(FileName,MovieNumbers(i),PatternNumber,SD(1));
    
    FigureProperties=struct('FigureNumber',3,'Subplot',[2 3 3],'TimeRange',[1 60],'AmplitudeRange',[-150 60],'FontSize',14,'Colors',['k' 'b' 'r' 'k' 'g' 'm' 'c'],'LineWidth',1,'YLabel','Signal [\muV]');
    y=NS512_PlotClustersOfSignaturesOnArrayLayoutWithMarksRotated(Traces,ChannelsPlot,WaveformTypes,ArrayID,FigureProperties,NS_GlobalConstants,[],[]);
end
%break;
for i = 1:length(MovieNumbers)
    Amplitudes =FindAmplitudes_PR( DataPath,1,MovieNumbers(i),[1:512],NS_GlobalConstants );
    CurrentAmplitude = Amplitudes(max(size(Amplitudes)));
    
    [DataTraces1,ArtifactDataTraces,DataChannels]=NS_ReadPreprocessedData(DataPath,DataPath,0,PatternNumber,MovieNumbers(i),0,0); 
    DataTraces0=DataTraces1(1:100,:,:)+350;
    SD=size(DataTraces0);
    ChannelTraces = DataTraces0(:,ChannelNumber,:);
    ChannelTraces2D = reshape(ChannelTraces,SD(1),SD(3));
    %subplot(2,2,i);
    h= plot(time,ChannelTraces2D'/0.27,x,y*Skalowanie*CurrentAmplitude-1000);
    set(h(101),'Color','Blue');
    set(h(101),'LineWidth',2);    

    WaveformTypes=NS_ReadClusterFile(FileName,MovieNumbers(i),PatternNumber,SD(1));
    artifactsIndex = find(WaveformTypes==1);
    spikeIndex = find(WaveformTypes==2);
    exclude=find(WaveformTypes>2);
    set(h(exclude),'Visible','off');
    set(h(artifactsIndex),'Color','Black');
    set(h(spikeIndex),'Color','Red');
    text(1.53,-820,['Amplitude = ' strcat(sprintf('%0.2g',CurrentAmplitude),' \muA')],'Fontsize',20);
    text(1.53,-950,['Efficacy = ' strcat(sprintf('%0.2g',length(spikeIndex)),' %')],'Fontsize',20);
    h1=gca;
    set(h1,'XLim',[-0.35 2.5]);
    set(h1,'YLim',[-1200 600]);

    h2=gca;
    set(h2,'Fontsize',20);
    set(h2,'XTick',[0:0.5:3]);
    set(h2,'XTickLabel',{'0' '0.5' '1' '1.5' '2' '2.5' '3'});
    set(h2,'YTick',[-1000:250:500]);
    set(h2,'YTickLabel',{'-1000' '-750' '-500' '-250' '0' '250' '500'});

    grid on
    %if i==3
        ylabel('Signal [\muV]','FontSize',20);
        xlabel('Time [ms]','FontSize',20);
    %else
    %    ylabel('')
    %end
    PlotName=['EffVsStim_amp' num2str(i) '_movie' num2str(MovieNumbers(i))];
    FullName = ['D:\Home\Rydygier\NEURO\analysis\retina\2009-11-27-0\data001\April2010\figures\ToMEA2010\' PlotName]

    hc = figure(1);
    set(hc,'PaperUnits','inches');
    set(hc,'PaperSize',[10 7]);
    set(hc,'PaperPosition',[0 0 10 7]);    
    print(hc, '-dtiff', '-r120', FullName);
end
break;
FullName = 'D:\Home\Rydygier\NEURO\analysis\retina\2009-11-27-0\data001\April2010\figures\toPaper\EfficVsStimCurr_font20';

hc = figure(1);
set(hc,'PaperUnits','inches');
        set(hc,'PaperSize',[10 10]);
        set(hc,'PaperPosition',[0 0 10 10]);
        
        %print(hc, '-dtiff', '-r120', FullName);

       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Wydajnosc stymulacji dla impulsow 50us i 100us - porownanie

FileName1 = [DataPath filesep 'ClusterFile_001_el242_50us'];
k = 1;
for i=1:2:MaxMovieNumber
    WaveformTypes50us=NS_ReadClusterFile(FileName1,i,PatternNumber,100)
    Spikes50(k,1) = i;
    Spikes50(k,2) = sum(WaveformTypes50us==2)+sum(WaveformTypes50us==4)
    k = k+1;
end

FileName2 = [DataPath filesep 'ClusterFile_001_el242_100us'];
k = 1;
for i=2:2:MaxMovieNumber+1
    WaveformTypes100us=NS_ReadClusterFile(FileName2,i,PatternNumber,SD(1));
    Spikes100(k,1) = i;
    Spikes100(k,2) = sum(WaveformTypes100us==2)+sum(WaveformTypes100us==4);
    k = k+1;
end

Effic50 = Spikes50(:,2);
Effic100 = Spikes100(:,2);
AllMovies50 = Spikes50(:,1);
AllMovies100 = Spikes100(:,1);

Amplitudes50 =FindAmplitudes_PR2( DataPath,1,[1:2:151],[1:512],NS_GlobalConstants );
Amplitudes100 =FindAmplitudes_PR2( DataPath,1,[2:2:152],[1:512],NS_GlobalConstants );

f50=find(diff(Amplitudes50)>0)
f100=find(diff(Amplitudes100)>0.002)
cf_50=SigmoidalFit_PR(Amplitudes50,Effic50);
cf_100=SigmoidalFit_PR(Amplitudes100,Effic100);
%figure;
%plot(Effic50);

figure(21);
%clf;
subplot(1,2,1); hh=plot(cf_50);
set(hh,'LineWidth',2.5);
hold on;

hh=plot(Amplitudes50(f50), Spikes50(f50,2),'bd','Markersize',8);
set(gca,'Xscale','log');
    
    legend off;
    grid off;
axis([0.05 4 0 100]);
set(gca,'XTick',[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 2 3 4]);
set(gca,'XTickLabel',{'0.1','','0.3','','','','','','','1','','3',''});
set(gca,'YTick',[0 25 50 75 100]);
    set(gca,'YTickLabel',{'0','25','50','75','100'});
%    set(gca,'GridLineStyle','--');
xlabel(gca,'Amp. [\muA]','FontSize',20);
ylabel(gca,'Eff. [%]','FontSize',20);
set(gca,'FontSize',20);

subplot(1,2,2); hh=plot(cf_100);
set(hh,'LineWidth',2.5);
hold on;
hh=plot(Amplitudes100(f100), Spikes100(f100,2),'bd','Markersize',8);
set(gca,'Xscale','log');
    legend off;
    grid off;
axis([0.05 4 0 100]);
set(gca,'XTick',[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 2 3 4]);
set(gca,'XTickLabel',{'0.1','','0.3','','','','','','','1','','3',''});
set(gca,'YTick',[0 25 50 75 100]);
%set(gca,'YTickLabel',{'0','','50','','100'});
xlabel(gca,'Amp. [\muA]','FontSize',20);
%set(gca,'GridLineStyle','--');
ylabel(gca,'','Fontsize',1);
set(gca,'FontSize',20);
hold off;

FullName3 = 'D:\Home\Rydygier\NEURO\analysis\retina\2009-11-27-0\data001\April2010\figures\toPaper\Effic50and100us4';
hc = figure(21);
set(hc,'PaperUnits','inches');
        set(hc,'PaperSize',[10 4]);
        set(hc,'PaperPosition',[0 0 10 4]);
        print(hc, '-dtiff', '-r120', FullName3);
%break;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

    %h=plot(Amplitudes, Spikes(:,2),'bd','Markersize',6);
    %set(h,'LineWidth',0.1);
    %hold on;
    h=plot(cf_);
    hold on;
    
    set(gca,'Xscale','log');
    set(h,'LineWidth',2);
    legend off;
    
    h=plot(Amplitudes, Spikes(:,2),'bd','Markersize',6);
    
    axis([0.1 4 0 100]);
    set(gca,'XTick',[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 2 3 4]);
    %set(gca,'XTickLabel',{'0.1','','0.3','','','','','','','1','','3',''});
    
    set(gca,'YTick',[0 25 50 75 100]);
    set(gca,'YTickLabel',{'0','','50','','100'});
    set(gca,'GridLineStyle','-');
    grid on;
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

    hold off;
end

FullName2 = 'D:\Home\Rydygier\NEURO\analysis\retina\2009-11-27-0\data001\April2010\figures\toPaper\sigma_fit_vertical_log3';

hc = figure(36);
set(hc,'PaperUnits','inches');
        set(hc,'PaperSize',[2 9]);
        set(hc,'PaperPosition',[0 0 2 9]);
        
        %print(hc, '-dtiff', '-r120', FullName2);

% for i=1:64
%     [StimChannels,Amplitudes]=NS_StimulatedChannels(DataPath,i,1,[1:512],NS_GlobalConstants);
%     if find(StimChannels==299)
%         pattern=i
%     end
% end