NS_GlobalConstants=NS_GenerateGlobalConstants(500);

load D:\Home\Pawel\analysis\slices\SlicesTTX\2010-09-14-0\analysis_2013_10_02\files\EI1;
load D:\Home\Pawel\analysis\slices\SlicesTTX\2010-09-14-0\analysis_2013_10_02\files\dane; %MovieNumber,RepetitionNumber,PatternNumber,Latency,PatternID

MovieFilePath='D:\Home\Data\slices\2010-09-14-0\movie002';
%FullName='D:\Home\Pawel\analysis\slices\SlicesTTX\2010-09-14-0\analysis_2012_04_10\histograms5\GaussParameters.bin';
ImagePath='C:\home\Pawel\nauka\512stim_paper\FiguresProby\';
NeuronFilePath='C:\home\Pawel\nauka\analiza\SlicesTTX\2010-09-14-0\analysis_2012_04_10\Vision_output\2010-09-14-0\data002\data002.neurons';% - pokoj 023    
PreProcessedDataPath='J:\analysis\2010-09-14-0\data002_preproc';

GoodSpikes=find(dane(1,:)>0); % for some reason, the array "dane" generated originally by the NS512_SpikeTimesToStimulationPatterns_v2, has all zeros for some spikes.
Indexes=zeros(1,length(dane));
IndexesForGood=NS512_SpikeNumberInStimulatedBurst_v2(dane(:,GoodSpikes));
Indexes(GoodSpikes)=IndexesForGood;

Patterns=[339 371 409 395 463 491 32 64 47 30 45 356];
%Patterns=[339 63 409 395 463 491 32 64 47 30 45 356];

ExamplaryPattern=356;
movieIDList = NS512_GetMoviesWithPatternIDNew(MovieFilePath,ExamplaryPattern);
AmplitudesForPattern=zeros(1,length(movieIDList));
for i=1:length(movieIDList)
    AmplitudesForPattern(i)=NS_AmplitudesForPattern_512_1el(PreProcessedDataPath,[1:512],ExamplaryPattern,movieIDList(i),NS_GlobalConstants);
end

figure(11)
clf

LC=[0 0.33 0.66]+0.05;
BC=[0 0.17 0.66 0.83]+0.03;
LeftCoordinates=[LC(1) LC(2) LC(3) LC(1) LC(2) LC(3) LC(1) LC(2) LC(3) LC(1) LC(2) LC(3)];
BottomCoordinates=[BC(4) BC(4) BC(4) BC(3) BC(3) BC(3) BC(2) BC(2) BC(2) BC(1) BC(1) BC(1)];
SubplotWidth=0.28;
SubplotHeight=0.14;

for p=1:length(Patterns)
    WholeAreaCoordinates(1)=LeftCoordinates(p);
    WholeAreaCoordinates(2)=BottomCoordinates(p);
    WholeAreaCoordinates(3)=SubplotWidth;
    WholeAreaCoordinates(4)=SubplotHeight;
    
    %subplot('position',[0.05 0.05 0.56 0.92])
    Subsubplot(WholeAreaCoordinates,[0.01 0.01 0.98 0.42]);
    hold on;
    movieIDList = NS512_GetMoviesWithPatternIDNew(MovieFilePath,Patterns(p));
    NumberOfSpikes=zeros(1,length(movieIDList));
    AllSpikeTimes=[];
    for m=1:length(movieIDList)
        SpikesIDs=find(dane(1,:)==movieIDList(m) & dane(3,:)==Patterns(p));
        %gdfgdfh=dane(:,SpikesIDs)
        %Indexes=NS512_SpikeNumberInStimulatedBurst(dane(:,SpikesIDs))
        SpikeTimes=dane(4,SpikesIDs);  
        SpikeIndexes=Indexes(SpikesIDs);
        AllSpikeTimes=[AllSpikeTimes SpikeTimes];
        NumberOfSpikes(m)=length(SpikeTimes);                        
        for sp=1:length(SpikeTimes)
            h1=plot([SpikeTimes(sp)/20 SpikeTimes(sp)/20],-[-m -m-0.5],'b-');
            SpikeNumberInBurst=SpikeIndexes(sp);
            if SpikeNumberInBurst==2
                set(h1,'Color','r');
            end
            set(h1,'LineWidth',1);
        end
        h1=gca;
        set(h1,'XLim',[0 30]);
        set(h1,'XTick',[0:5:30]);
        MoviesToMark=[2:5:17];
        set(h1,'YTick',MoviesToMark);
        %for m3=1:length(MoviesToMark)
        %    MovietoMark=MoviesToMark(m3);
        %    f1=get(h1);
            
            
        set(h1,'YTick',[2:5:17]);
        set(h1,'YTickLabel',{'0.33' '0.53' '0.86' '1.41'});
        if p==10
            xlabel('Time [ms]');
            ylabel('Amp [\muA]');
            %set(h1,'YTickLabel',{'0.4' '0.7' '1.2' '2.3'});
        %else
        %    set(h1,'XTickLabel','');
        %    set(h1,'YTickLabel','');
        end
    end    
    Subsubplot(WholeAreaCoordinates,[0.01,0.56,0.46,0.42]);
    plot(int32(NumberOfSpikes/102*100),'bd-')
    axis([1 17 0 125]);
    grid on
    h3=text(2,100,num2str(p));
    set(h3,'Color','r');
    h4=gca;
    set(h4,'XTick',[1:2:17]);
    set(h4,'YTick',[0:25:125]);
    %if p==10
    %set(h4,'XTickLabel',{ '' '0.4' '' '0.7' '' '1.2' '' '2.3' ''});
    set(h4,'XTickLabel',{'' '0.36' '' '0.53' '' '0.78' '' '1.10' ''});
    set(h4,'XTickLabel',{'' '0.36' '' '' '0.63' '' '' '1.10' ''});
    %set(h4,'XTickLabel',{'0.30' '' '0.43' '' '0.63' '' '0.93' '' '1.41'});
    set(h4,'YTickLabel',{'0' '' '50' '' '100' ''});
    
    if p==10
        xlabel('Amp [\muA]');
        ylabel('Eff [%]');
    end       
    
    Subsubplot(WholeAreaCoordinates,[0.53,0.56,0.46,0.42])
    LatenciesForGivenElectrode=AllSpikeTimes;
    p1=hist(LatenciesForGivenElectrode/20,[5:10:600]/20);
    hist(LatenciesForGivenElectrode/20,[5:10:600]/20);
    
    %[CFs,leg,moje]=NS512_FitWithMultiGauss_2015_05_26([1:60],p1);
    CFs=NS512_FitWithMultiGauss([1:60],p1);
    hold on
    Peak=zeros(60,1);
    for peak=1:length(CFs)
        FL=CFs{peak};
        Peak=Peak+FL([1:60]);        
    end
        h=plot([5:10:600]/20,Peak,'r-');
        grid on;
        set(h,'LineWidth',2)
        h=gca;
        set(h,'XLim',[0 600]/20);
        set(h,'YLim',[0 max(p1)*1.5]);
        set(h,'YTick',[0 0.5 1 1.5]*max(p1));
        set(h,'XTick',[0:5:30]);
        set(h,'YTickLabel','');
        set(h,'XTickLabel',{'0' '' '10' '' '20' '' '30'});
        if p==10
            %set(h,'XTickLabel',{'0' '' '10' '' '20' '' '30'});
            xlabel('Time [ms]');
        else
            %set(h,'XTickLabel','');
        end
        %ylabel('N');
        %text(450,40,['A=' num2str(CFs{1}.A,'%8.2f')]);
        text(17,max(p1)*1.3,['\tau = ' num2str(CFs{1}.tau*10/20,'%8.1f')]);
        text(17,max(p1)*1,['\sigma = ' num2str(CFs{1}.sigma*10/20,'%8.1f')]);                
    %end
    if p==11        
        figure(101);
        plot(int32(NumberOfSpikes/102*100),'bd-');
        h=gca;
        set(h,'FontSize',24);
        h=gcf;
        %FullImageName=['C:\home\Pawel\nauka\analiza\SlicesTTX\2010-09-14-0\analysis_2013_10_02\figures\figure1_eff.tif'];
        FullImageName=[ImagePath 'figure2_eff.tif'];
        set(h,'PaperUnits','inches');
        set(h,'PaperSize',[8 5]*2);
        set(h,'PaperPosition',[0 0 8 5]*2); 
        print(h, '-dtiff', '-r120', FullImageName);
        
        figure(102);
        hist(LatenciesForGivenElectrode,[5:10:600]);
        hold on;
        h=plot([5:10:600],Peak,'r-');
        set(h,'LineWidth',2)
        h=gca;
        set(h,'XLim',[0 600]);
        set(h,'FontSize',24);
        xlabel('Time [ms]');
        ylabel('N');
        %text(450,40,['A=' num2str(CFs{1}.A,'%8.2f')]);
        text(20,35,['\tau = ' num2str(CFs{1}.tau*10/20,'%8.1f') ' ms']);
        text(20,30,['\sigma = ' num2str(CFs{1}.sigma*10/20,'%8.1f') ' ms']);
        h=gcf
        FullImageName=[ImagePath 'figure2_hist.tif'];
        set(h,'PaperUnits','inches');
        set(h,'PaperSize',[8 5]*2);
        set(h,'PaperPosition',[0 0 8 5]*2); 
        print(h, '-dtiff', '-r120', FullImageName);
        
        figure(103);
        clf;
        hold on;
                
        for m=1:length(movieIDList)
            m
            SpikesIDs=find(dane(1,:)==movieIDList(m) & dane(3,:)==Patterns(p));       
            SpikeTimes=dane(4,SpikesIDs) 
            SpikeIndexes=Indexes(SpikesIDs);
            AllSpikeTimes=[AllSpikeTimes SpikeTimes];
            NumberOfSpikes(m)=length(SpikeTimes);                        
            for sp=1:length(SpikeTimes)
                h1=plot([SpikeTimes(sp)/20 SpikeTimes(sp)/20],-[-m -m-0.5],'b-');
                SpikeNumberInBurst=SpikeIndexes(sp);
                if SpikeNumberInBurst>=2
                    set(h1,'Color','r');
                end
                set(h1,'LineWidth',1);
            end
            
            AmplitudesToLabelIndexes=[1:2:17];
            AmplitudesToLabel=[1:2:17];
            for a1=1:length(AmplitudesToLabelIndexes)
                [Amplitudes]=NS_AmplitudesForPattern_512('J:\analysis\2010-09-14-0\data002_preproc',[1:512],13,AmplitudesToLabelIndexes(a1)*8,NS_GlobalConstants);
                AmplitudesToLabel(a1)=max(abs(Amplitudes));
                bleble{a1}=num2str(max(abs(Amplitudes)),'%8.2f');
            end
                
            h1=gca;
            set(h1,'XLim',[0 30]);
            set(h1,'XTick',[0:5:30]);
            set(h1,'YTick',[AmplitudesToLabelIndexes]);
            set(h1,'YTickLabel',bleble)
            %set(h1,'YTickLabel',num2str(AmplitudesToLabel));
            %if p==10
                %xlabel('Time [ms]');
                %ylabel('Amp [\muA]');
            %set(h1,'YTickLabel',{'0.4' '0.7' '1.2' '2.3'});
        %else
        %    set(h1,'XTickLabel','');
        %    set(h1,'YTickLabel','');
            %end
        end   
        h2=xlabel('Time [ms]');
        set(h2,'FontSize',24);
        h2=ylabel('Amplitude [\muA]');
        set(h2,'FontSize',24)
        h2=gca
        set(h2,'FontSize',24)
        
        h=gcf;
        FullImageName=[ImagePath 'figure2_raster.tif'];
        set(h,'PaperUnits','inches');
        set(h,'PaperSize',[8 5]*2);
        set(h,'PaperPosition',[0 0 8 5]*2); 
        print(h, '-dtiff', '-r120', FullImageName);
                
        figure(11)
    end        
end

figure(11)
subplot('position',[0.15 0.39 0.7 0.22])
h=NS512_ShowEIAsCircles4(EI1/1.5,500,[1:512],[],Patterns,[-1000 1000],[-500 500]);
set(h,'Visible','off');

FullImageName=[ImagePath 'figure5.tif'];
h=gcf;
set(h,'PaperUnits','inches');
set(h,'PaperSize',[10 12]);
set(h,'PaperPosition',[0 0 10 12]); 
print(h, '-dtiff', '-r240', FullImageName);

figure(104)
clf
h=NS512_ShowEIAsCircles5(EI1/1.5,500,[1:512],[],Patterns(11),[-1000 1000],[-500 500]);
set(h,'Visible','off');
FullImageName=[ImagePath 'figure1_array.tif'];
h=gcf;
set(h,'PaperUnits','inches');
set(h,'PaperSize',[10 5]);
set(h,'PaperPosition',[0 0 10 5]); 
print(h, '-dtiff', '-r240', FullImageName);