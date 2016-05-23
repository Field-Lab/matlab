%clear
NS_GlobalConstants=NS_GenerateGlobalConstants(500);
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(500);

MovieFile='D:\Home\Pawel\analysis\slices\2013\2013-12-12-3-PH\MovieFiles\movie001';
DataPath='D:\Home\Pawel\analysis\slices\2013\2013-12-12-3-PH\data001preproc'
ArtifactDataPath=DataPath;
FigurePath='D:\Home\Pawel\analysis\slices\2013\2013-12-12-3-PH\data001_Matlab\SpikeAnalysis\2012-12-12-3-PH_001';

figure(2);
clf;

NumberOfStimulatedElectrodes=zeros(NumberOfAmplitudes,512);
StimulatedElectrodes=zeros(NumberOfAmplitudes,512,512);
AllPatterns=[]
for Amplitude=2:NumberOfAmplitudes
    Amplitude
    tic
    for MovieSeqeunce=1:NumberOfMovieSequences
        Movie=(Amplitude-1)*NumberOfMovieSequences2+MovieSeqeunce
        ChunkData=NS_MovieData_GlobalPath(MovieFile,Movie,NS_GlobalConstants);
        [PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,MovieData]=NS_DecodeMovieDataChunk(ChunkData);  
        
        PatternsForMovie=unique(MovieData(2:3:length(MovieData)));
        break
        AllPatterns=[AllPatterns' PatternsForMovie']';
        for p=1:length(PatternsForMovie)
            Pattern=PatternsForMovie(p)            
            [DataTraces,ArtifactDataTraces,Channels]=NS_ReadPreprocessedData(DataPath,ArtifactDataPath,0,Pattern,Movie,1500,0);
            
            Events=NS512_DetectSpikes(DataTraces,40,10);
            SpikesPerRecordingElectrode=sum(Events);
            [a,b]=find(Events>0);
            n=hist(b,[1:512]);
            ElectrodesWithSpikes=find(n>25);
            NumberOfStimulatedElectrodes(Amplitude,Pattern)=length(ElectrodesWithSpikes);
            StimulatedElectrodes(Amplitude,Pattern,:)=n;
            %break;
            ElectrodesWithSpikes2=find(n>10);
            clear X;
            clear Y;
            clear X2;
            clear Y2;
            for i=1:length(ElectrodesWithSpikes)
                X(i)=electrodeMap.getXPosition(ElectrodesWithSpikes(i));
                Y(i)=electrodeMap.getYPosition(ElectrodesWithSpikes(i));        
            end
            
            %for i=1:length(ElectrodesWithSpikes2)
            %    X2(i)=electrodeMap.getXPosition(ElectrodesWithSpikes2(i));
            %    Y2(i)=electrodeMap.getYPosition(ElectrodesWithSpikes2(i));        
            %end
            
            PatternX=electrodeMap.getXPosition(Pattern);
            PatternY=electrodeMap.getYPosition(Pattern);
            
            SubplotColumn=(PatternX+945)/120;
            SubplotRow=(PatternY+450)/120;
            
            %figure(1)
            subplot('position',[SubplotColumn*0.06+0.02 SubplotRow*0.06+0.52 0.05 0.05]);   
            hold on
            %if ElectrodesWithSpikes2
            %    h1=plot(X2,Y2,'bo');
            %    set(h1,'MarkerFaceColor','b');
            %   set(h1,'MarkerSize',2)
            %end
                       
            if ElectrodesWithSpikes
                h1=plot(X,Y,'ro');
                set(h1,'MarkerFaceColor','r');
                set(h1,'MarkerSize',4)
            end
            
            h1=plot(PatternX,PatternY,'gd');    
            set(h1,'MarkerFaceColor','g');
            axis([-1000 1000 -500 500]);
            h=gca;
            set(h,'XTickLabel','');
            set(h,'YTickLabel',''); 
            set(h,'Box','on');
            h=text(450,-360,num2str(Pattern));
            set(h,'FontSize',8);
            
            %figure(2)
            subplot('position',[SubplotColumn*0.06+0.02 SubplotRow*0.06+0.02 0.05 0.05]);    
            h=gca;
            set(h,'XTickLabel','');
            set(h,'YTickLabel','');                                      
            %hold on;
            if 0==1%ElectrodesWithSpikes2
                if ElectrodesWithSpikes
                    h=plot([1:50],sum(Events(:,ElectrodesWithSpikes2)/length(ElectrodesWithSpikes2),2),'bd-',[1:50],sum(Events(:,ElectrodesWithSpikes)/length(ElectrodesWithSpikes),2),'rd-');
                    set(h(1),'MarkerFaceColor','b');
                    set(h(2),'MarkerFaceColor','r');
                else
                    h=plot([1:50],sum(Events(:,ElectrodesWithSpikes2)/length(ElectrodesWithSpikes2),2),'bd-');
                    set(h,'MarkerFaceColor','b');
                end                
                set(h,'MarkerSize',3);
                axis([0 50 0 1]);
                h=gca
                set(h,'XTickLabel','')
                set(h,'YTickLabel','')  
                set(h,'Box','on');
                %grid on;
                set(h,'XTick',[0:5:50]);
                set(h,'YTick',[0:0.2:1]);
                h=text(40,0.1,num2str(Pattern));
                set(h,'FontSize',8);
                %grid on;
            end                                         
        end       
    end
    FullName=[FigurePath '\Amplitude' num2str(Amplitude) '.tif']; 
    h=gcf;
    set(h,'PaperUnits','inches');
    set(h,'PaperSize',[16 9]);
    set(h,'PaperPosition',[0 0 16 9]); 
    %print(h, '-dtiff', '-r120', FullName);
    toc
end