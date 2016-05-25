% Detekcja spikow z progiem 40 (czyli 150 uV). 

electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(500); %define the electrode map - must be different than 1 for silicon probes
NS_GlobalConstants=NS_GenerateGlobalConstants(500);

% 1) Spike detection
DataSet=1;
%Stim512PaperFinal_DetectAllSpikes3_rat.m

DataSet=2;
%Stim512PaperFinal_DetectAllSpikes3_rat.m

%TimeWindow=10;
%DataSet=1;
%Stim512PaperFinal_fig4_obliczenia;
%DataSet=2;
%Stim512PaperFinal_fig4_obliczenia;

Radius=1;

X=zeros(1,512);
Y=X;
for i=1:512
    X(i)=electrodeMap.getXPosition(i);
    Y(i)=electrodeMap.getYPosition(i);
end

figure(2)
clf
%hold on
N=0;

NumberOfDataSets=2;
FontSize=18;
LineWidth=2;
FontSize2=24;

for DataSet=1:2
    switch DataSet
        case 1 % without drugs
            %GaussesFilesPath='D:\Home\Pawel\analysis\slices\512stim_paper\SpikesAnalysis\2010-09-14-0\GaussesFiles2\';
            %GaussesFilesPath='C:\home\Pawel\nauka\512stim_paper\SpikesAnalysis\2010-09-14-0\data002\GaussesFiles-2015-05-24\';
            GaussesFilesPath='C:\home\Pawel\nauka\512stim_paper\SpikesAnalysis\2010-09-14-0\data002\GaussesFiles-2015-06-11\';
            MovieFilePath='C:\home\Pawel\nauka\512stim_paper\MovieFiles\2010-09-14-0\movie002';
        case 2 % with drugs
            %GaussesFilesPath='C:\home\Pawel\nauka\512stim_paper\SpikesAnalysis\2010-09-14-0\data005\GaussesFiles-2015-05-24\';   
            GaussesFilesPath='C:\home\Pawel\nauka\512stim_paper\SpikesAnalysis\2010-09-14-0\data005\GaussesFiles-2015-06-11\';
            MovieFilePath='C:\home\Pawel\nauka\512stim_paper\MovieFiles\2010-09-14-0\movie005';             
        otherwise
            error('glaglaglabulbulbul!!');
    end
    [NumberOfMovies,NumberOfPatternsPerMovie,AllPatternsUsed,Patterns]=NS512_MoviePatterns(MovieFilePath,NS_GlobalConstants);
    
    h1=[];
    for Pattern=AllPatternsUsed%(1:10)
        fid=fopen([GaussesFilesPath 'AllGausses_p' num2str(Pattern)],'r');
        a=fread(fid,'double');
        fclose(fid);
    
        b0=reshape(a,length(a)/6,6);  
        GoodFit=find(b0(:,6)<=0.1);
        b=b0(GoodFit,:);
        sb=size(b)
        ChannelsToExclude=electrodeMap.getAdjacentsTo(Pattern,Radius)';
        %g=zeros(1,1);
        g=[];
        for c=ChannelsToExclude
            l=find(b(:,2)==c);
            g=[g' l']';
        end            
        HistogramsToPlot=setdiff([1:sb(1)],g);            
        LowJitter=find(b(HistogramsToPlot,5)<9);
    
        Distances=sqrt((X(b(HistogramsToPlot,1))- X(b(HistogramsToPlot,2))).^2+(Y(b(HistogramsToPlot,1))- Y(b(HistogramsToPlot,2))).^2);    N=N+length(HistogramsToPlot);
        %HistogramsToExclude=find(b(:,2)==ChannelsToExclude)
        %subplot(NumberOfDataSets,4,(DataSet-1)*4+2)
        subplot('Position', [0.32, 0.58-(DataSet-1)*0.51, 0.15, 0.32]);
        hold on
        plot(b(HistogramsToPlot,4)/20,b(HistogramsToPlot,5)/20,'bd');
        plot(b(HistogramsToPlot(LowJitter),4)/20,b(HistogramsToPlot(LowJitter),5)/20,'rd');
        axis([0 30 0 6])
        h=gca;
                set(h,'LineWidth',LineWidth);
        set(h,'Fontsize',FontSize);
        %if DataSet==NumberOfDataSets
            xlabel('Latency [ms]');
            ylabel('Jitter [ms]');
        %end
        grid on
        h11=gca;
        set(h11,'Fontsize',FontSize);
        ax2 = axes('Position', [0.28, 0.95-(DataSet-1)*0.51, 0.04, 0.04]);
        axis([-1 1 -1 1])
        %plot([0],[0],'w');
        h5=gca;
        set(h5,'Box','off')
        set(h5,'Visible','off')
        switch DataSet
            case 1
                h6=text(0,0,'B)');
            case 2
                h6=text(0,0,'F)');
        end
        set(h6,'FontSize',FontSize2)
        
        %subplot(NumberOfDataSets,4,(DataSet-1)*4+3)
        subplot('Position', [0.57, 0.58-(DataSet-1)*0.51, 0.15, 0.32]);
        hold on
        h2=plot(Distances,b(HistogramsToPlot,4)/20,'bd');
        %set(h2,'MarkerSize',4);
        %set(h2,'MarkerFaceColor','b');
        plot(Distances(LowJitter),b(HistogramsToPlot(LowJitter),4)/20,'rd');
        axis([0 2000 0 30])
        h=gca;
                set(h,'LineWidth',LineWidth);
        set(h,'XTick',[0:500:2000]);
        set(h,'XTickLabel',{'0' '' '1000' '' '2000'})
        set(h,'Fontsize',FontSize);
        %if DataSet==NumberOfDataSets
            xlabel('Distance [\mum]');
            ylabel('Latency [ms]');
        %end
        grid on
        ax2 = axes('Position', [0.53, 0.95-(DataSet-1)*0.51, 0.04, 0.04]);
        axis([-1 1 -1 1])
        %plot([0],[0],'w');
        h5=gca;
        set(h5,'Box','off');
        set(h5,'Visible','off');
        switch DataSet
            case 1
                h6=text(0,0,'C)');
            case 2
                h6=text(0,0,'G)');
        end
        set(h6,'FontSize',FontSize2);
        
        %subplot(NumberOfDataSets,4,(DataSet-1)*4+4)
        subplot('Position', [0.82, 0.58-(DataSet-1)*0.51, 0.15, 0.32]);
        hold on
        plot(Distances,b(HistogramsToPlot,5)/20,'bd');
        plot(Distances(LowJitter),b(HistogramsToPlot(LowJitter),5)/20,'rd')
        axis([0 2000 0 6]);
        h=gca;
        set(h,'XTick',[0:500:2000]);
        set(h,'XTickLabel',{'0' '' '1000' '' '2000'});
        set(h,'Fontsize',FontSize);
        set(h,'LineWidth',LineWidth);
        %if DataSet==NumberOfDataSets
            xlabel('Distance [\mum]');
            ylabel('Jitter [ms]');
        %end
        grid on
        ax2 = axes('Position', [0.78, 0.95-(DataSet-1)*0.51, 0.04, 0.04]);
        axis([-1 1 -1 1])
        %plot([0],[0],'w');
        h5=gca;
        set(h5,'Box','off');
        set(h5,'Visible','off');
        switch DataSet
            case 1
                h6=text(0,0,'D)');
            case 2
                h6=text(0,0,'H)');
        end
        set(h6,'FontSize',FontSize2);
        
        h1=[h1 b(HistogramsToPlot,5)'];
    end
    %subplot(NumberOfDataSets,4,(DataSet-1)*4+1);
    subplot('Position', [0.06, 0.58-(DataSet-1)*0.51, 0.15, 0.32]);
    hist(h1/20,[1:600]/20);
    h=gca;
    set(h,'XLim',[0 6]);
    set(h,'Fontsize',FontSize);
    set(h,'Box','off');
    set(h,'LineWidth',LineWidth);
    h12=xlabel('Jitter [ms]');
    set(h12,'Fontsize',FontSize);
    h12=ylabel('Counts');
    set(h12,'Fontsize',FontSize);
    %grid on
    ax2 = axes('Position', [0.13, 0.68-(DataSet-1)*0.51, 0.09, 0.24]);
    hist(ax2,h1/20,[1:600]/20)
    axis([0 3 0 200])
    grid on
    h=gca;
    set(h,'LineWidth',LineWidth);
    set(h,'Fontsize',FontSize);
    h12=xlabel('Jitter [ms]');
    set(h12,'Fontsize',FontSize);
    h12=ylabel('Counts');
    set(h12,'Fontsize',FontSize);
    %subplot('position',[0.2, 0.7, 0.1, 0.1])
    ax2 = axes('Position', [0.03, 0.95-(DataSet-1)*0.51, 0.04, 0.04]);
    axis([-1 1 -1 1]);
    %plot([0],[0],'w');
    h5=gca;
    set(h5,'Box','off');
    set(h5,'Visible','off');
    switch DataSet
        case 1
            h6=text(0,0,'A)');
        case 2
            h6=text(0,0,'E)');
    end
    set(h6,'FontSize',FontSize2)
end

h=gcf;
set(h,'PaperUnits','inches');
set(h,'PaperSize',[16 9]);
set(h,'PaperPosition',[0 0 16 9]); 
FullFigurePath='C:\home\Pawel\nauka\512stim_paper\FiguresProby\figure4_7.tif';
print(h, '-dtiff', '-r120', FullFigurePath); 