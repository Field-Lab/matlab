%clear
NS_GlobalConstants=NS_GenerateGlobalConstants(500);
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(500);

MovieFile='D:\Home\Pawel\analysis\slices\2013\2013-12-12-3-PH\MovieFiles\movie001';
DataPath='D:\Home\Pawel\analysis\slices\2013\2013-12-12-3-PH\data001preproc'
ArtifactDataPath=DataPath;
%FigurePath='D:\Home\Pawel\analysis\slices\2013\2013-12-12-3-PH\data001_Matlab\SpikeAnalysis\2012-12-12-3-PH_001';

figure(2);
clf;

electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(500); %define the electrode map - must be different than 1 for silicon probes

for i=1:512
    X(i)=electrodeMap.getXPosition(i);
    Y(i)=electrodeMap.getYPosition(i);        
end

%MapaKolorow=jet;

NumberOfStimulatedElectrodes=zeros(NumberOfAmplitudes,512);
StimulatedElectrodes=zeros(NumberOfAmplitudes,512,512);
AllPatterns=[]
for Amplitude=18%2:NumberOfAmplitudes
    Amplitude
    tic
    for MovieSeqeunce=1:NumberOfMovieSequences
        Movie=(Amplitude-1)*NumberOfMovieSequences2+MovieSeqeunce
        ChunkData=NS_MovieData_GlobalPath(MovieFile,Movie,NS_GlobalConstants);
        [PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,MovieData]=NS_DecodeMovieDataChunk(ChunkData);  
        
        PatternsForMovie=unique(MovieData(2:3:length(MovieData)));
        AllPatterns=[AllPatterns' PatternsForMovie']';
        for p=1:length(PatternsForMovie)
            Pattern=PatternsForMovie(p)            
            [DataTraces,ArtifactDataTraces,Channels]=NS_ReadPreprocessedData(DataPath,ArtifactDataPath,0,Pattern,Movie,1500,0);
            
            Events=NS512_DetectSpikes(DataTraces,40,10);
            [a,b]=find(Events>0);
            n=hist(b,[1:512]);
            
            PatternX=electrodeMap.getXPosition(Pattern);
            PatternY=electrodeMap.getYPosition(Pattern);
            
            SubplotColumn=(PatternX+945)/120;
            SubplotRow=(PatternY+450)/120;
            
            ConvertedData=ChannelDataToColors2(n,X,Y);
            subplot('position',[SubplotColumn*0.06+0.02 SubplotRow*0.06+0.52 0.05 0.05]);   
            image(ConvertedData')
            colormap(jet(50))
            
            hold on
            
            h10=plot(X(Pattern)/30+33,17-round((Y(Pattern)+490)/60),'go');
            set(h10,'MarkerSize',8)            
            %set(h10,'MarkerFaceColor','g');
            h=gca;
            set(h,'XTickLabel','');
            set(h,'YTickLabel',''); 
            %et(h,'Box','on');
            %=text(450,-360,num2str(Pattern));
            %et(h,'FontSize',8);                                                                                                   
        end       
    end
    FullName=[FigurePath '\Amplitude' num2str(Amplitude) '.tif']; 
    h=gcf;
    set(h,'PaperUnits','inches');
    set(h,'PaperSize',[16 9]);
    set(h,'PaperPosition',[0 0 16 9]); 
    print(h, '-dtiff', '-r120', FullName);
    toc
end