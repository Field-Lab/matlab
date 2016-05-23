electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(500); %define the electrode map - must be different than 1 for silicon probes
NS_GlobalConstants=NS_GenerateGlobalConstants(500);

Radius=1;

X=zeros(1,512);
Y=X;
for i=1:512
    X(i)=electrodeMap.getXPosition(i);
    Y(i)=electrodeMap.getYPosition(i);
end

figure(1)
clf
%hold on
N=0;

NumberOfDataSets=2;

for DataSet=1  %:2
    switch DataSet
        case 1
            GaussesFilesPath='D:\Home\Pawel\analysis\slices\512stim_paper\SpikesAnalysis\2010-09-14-0\GaussesFiles2\';
            MovieFilePath='D:\Home\Pawel\analysis\slices\512stim_paper\SpikesAnalysis\2010-09-14-0\movie002';
        case 2
            GaussesFilesPath='D:\Home\Pawel\analysis\slices\512stim_paper\SpikesAnalysis\2010-09-14-0\data005\GaussesFiles\';
            %GaussesFilesPath='D:\Home\Pawel\analysis\slices\512stim_paper\SpikesAnalysis\2010-09-14-0\GaussesFiles\';
            %MovieFilePath='D:\Home\Pawel\analysis\slices\512stim_paper\SpikesAnalysis\2010-09-14-0\movie002';
            MovieFilePath='G:\analysis\slices\2010-09-14-0\TTX_sub\movie005';
        %case 1
        %    GaussesFilesPath='D:\Home\Pawel\analysis\slices\512stim_paper\SpikesAnalysis\2013-12-12-0\data001\GaussesFiles\';
        %   MovieFilePath='J:\data\2013-12-12-0-PH\movie001';
        %case 2
        %    GaussesFilesPath='D:\Home\Pawel\analysis\slices\512stim_paper\SpikesAnalysis\2013-12-12-0\data004\GaussesFiles\';
        %    MovieFilePath='J:\data\2013-12-12-0-PH\movie004';
        %case 1
            %GaussesFilesPath='D:\Home\Pawel\analysis\slices\512stim_paper\SpikesAnalysis\2013-12-12-3\data001\GaussesFiles\';
            %MovieFilePath='J:\data\2013-12-12-3-PH\movie001';
        %case 2
            %GaussesFilesPath='D:\Home\Pawel\analysis\slices\512stim_paper\SpikesAnalysis\2013-12-15-0\data004\GaussesFiles\';
            %MovieFilePath='J:\data\2013-12-12-0-PH\movie004';        
        otherwise
            error('glaglaglabulbulbul!!');
    end
    [NumberOfMovies,NumberOfPatternsPerMovie,AllPatternsUsed,Patterns]=NS512_MoviePatterns(MovieFilePath,NS_GlobalConstants);
    
    h1=[];
    for Pattern=AllPatternsUsed%(1:10)
        fid=fopen([GaussesFilesPath 'AllGausses_p' num2str(Pattern)],'r');
        a=fread(fid,'double');
        fclose(fid);
    
        b=reshape(a,length(a)/5,5);  
        sb=size(b)
        ChannelsToExclude=electrodeMap.getAdjacentsTo(Pattern,Radius)';
        %g=zeros(1,1);
        g=[];
        for c=ChannelsToExclude
            l=find(b(:,2)==c);
            g=[g' l']'
        end            
        HistogramsToPlot=setdiff([1:sb(1)],g);            
        LowJitter=find(b(HistogramsToPlot,5)<10)
    
        Distances=sqrt((X(b(HistogramsToPlot,1))- X(b(HistogramsToPlot,2))).^2+(Y(b(HistogramsToPlot,1))- Y(b(HistogramsToPlot,2))).^2);    N=N+length(HistogramsToPlot);
        %HistogramsToExclude=find(b(:,2)==ChannelsToExclude)
        subplot(NumberOfDataSets,4,(DataSet-1)*4+2)
        hold on
        plot(b(HistogramsToPlot,4)/20,b(HistogramsToPlot,5)/20,'bd');
        plot(b(HistogramsToPlot(LowJitter),4)/20,b(HistogramsToPlot(LowJitter),5)/20,'rd');
        axis([0 30 0 6])
        %if DataSet==NumberOfDataSets
            xlabel('Latency [ms]');
            ylabel('Jitter [ms]');
        %end
        grid on
        
        subplot(NumberOfDataSets,4,(DataSet-1)*4+3)
        hold on
        h2=plot(Distances,b(HistogramsToPlot,4)/20,'bd');
        %set(h2,'MarkerSize',4);
        %set(h2,'MarkerFaceColor','b');
        plot(Distances(LowJitter),b(HistogramsToPlot(LowJitter),4)/20,'rd');
        axis([0 2000 0 30])
        %if DataSet==NumberOfDataSets
            xlabel('Distance [\mum]');
            ylabel('Latency [ms]');
        %end
        grid on
        
        subplot(NumberOfDataSets,4,(DataSet-1)*4+4)
        hold on
        plot(Distances,b(HistogramsToPlot,5)/20,'bd')
        plot(Distances(LowJitter),b(HistogramsToPlot(LowJitter),5)/20,'rd')
        axis([0 2000 0 6])
        %if DataSet==NumberOfDataSets
            xlabel('Distance [\mum]');
            ylabel('Jitter [ms]');
        %end
        grid on
        h1=[h1 b(HistogramsToPlot,5)'];
    end
    subplot(NumberOfDataSets,4,(DataSet-1)*4+1);
    hist(h1/20,[1:600]/20)
    h=gca
    set(h,'XLim',[0 6])
    xlabel('Jitter [ms]');
    grid on
end

h=gcf;
set(h,'PaperUnits','inches');
set(h,'PaperSize',[16 9]);
set(h,'PaperPosition',[0 0 16 9]); 
FullFigurePath='C:\home\Pawel\nauka\512stim_paper\FiguresProby\figure4.tif';
print(h, '-dtiff', '-r100', FullFigurePath); 