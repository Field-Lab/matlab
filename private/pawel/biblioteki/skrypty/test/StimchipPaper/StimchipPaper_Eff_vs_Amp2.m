clear;
NS_GlobalConstants=NS_GenerateGlobalConstants(61);
DataPath= 'G:\komp_NS512\2010\analysis\retina_61\2010-09-21-0\data003\';

DataPath='E:\Stimchip_Paper_data\data003';
FileName = [DataPath filesep 'ClusterFile_003_StimPaper'];
FontSize=20;
FrameThickness=2;

ChannelNumber = 16;
Patterns=[16 13 14 1 60 56 55];
Movies=[1:26];

figure(2)
clf;
for i=1:length(Patterns)
    Pattern=Patterns(i)
    for j=1:length(Movies)
        Movie=Movies(j);
        WaveformTypes=NS_ReadClusterFile(FileName,j,Pattern,50);
        Efficacy(j)=sum(WaveformTypes==2)*2;
    end
    Amplitudes=FindAmplitudes_PR2(DataPath,1,Movies,[1:64],NS_GlobalConstants);
    SF=SigmoidalFit_PR(Amplitudes,Efficacy);
    figure(2)
    %subplot(length(Patterns),1,i)
    subplot('Position',[0.3 1.01-i*0.135 0.6 0.12*0.92]);
    h2=plot(SF,'r-');
    set(h2,'LineWidth',6);
    hold on;
    %h1=plot(Amplitudes,Efficacy,'bd');
    
    %h1=plot(SF,'r-',Amplitudes,Efficacy,'bd');
    %set(h1(1),'LineWidth',6);
    grid on;
    legend off;
    h=text(0.18,80,num2str(i));
    set(h,'FontSize',FontSize);
    set(h,'BackgroundColor','w')            
    
    h=gca;
    set(h,'FontSize',FontSize);
    set(h,'LineWidth',FrameThickness);
    set(h,'XScale','log');
    set(h,'YTick',[0:25:100]);
    set(h,'XLim',[0.15 1.62]);
    set(h,'XTick',[0.15 0.2:0.1:1.6]);
    if i==length(Patterns)
        set(h,'XTickLabel',{'0.15' '' '' '' '0.5' '' '' '' '' '' '' '' '' '' '1.5' ''});
        xlabel('Amplitude [\muA]');
    else        
        set(h,'XTickLabel',[]);
        h=xlabel('');
        set(h,'FontSize',FontSize);
    end
    h=ylabel('\epsilon [%]');
    set(h,'FontSize',FontSize);
    p1=get(h,'Position')
    p1(1)=0.07;
    set(h,'Position',p1);
    hold on;    
end
break
h=gcf;
FullName=['C:\home\pawel\nauka\Stimchip_paper\obrazki\Eff_vs_Amp'];
set(h,'PaperUnits','inches');
set(h,'PaperSize',[2.7 12]);
set(h,'PaperPosition',[0 0 2.7 12]);  
print(h, '-dtiff', '-r400', FullName);