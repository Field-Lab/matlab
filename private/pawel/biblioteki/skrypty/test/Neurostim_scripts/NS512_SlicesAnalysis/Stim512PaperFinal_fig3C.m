%clear
NS_GlobalConstants=NS_GenerateGlobalConstants(500);
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(500);

electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(500); %define the electrode map - must be different than 1 for silicon probes

for i=1:512
    X(i)=electrodeMap.getXPosition(i);
    Y(i)=electrodeMap.getYPosition(i);        
end

NumberOfStimulatedElectrodes=zeros(NumberOfAmplitudes,512);
AllPatterns=[]
for Amplitude=18  
    for MovieSeqeunce=1:NumberOfMovieSequences
        Movie=(Amplitude-1)*NumberOfMovieSequences2+MovieSeqeunce
        ChunkData=NS_MovieData_GlobalPath(MovieFile,Movie,NS_GlobalConstants);
        [PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,MovieData]=NS_DecodeMovieDataChunk(ChunkData);  
        
        PatternsForMovie=unique(MovieData(2:3:length(MovieData)));
        AllPatterns=[AllPatterns' PatternsForMovie']';
        for p=1:length(PatternsForMovie)
            Pattern=PatternsForMovie(p);            
            n=StimulatedElectrodesMouse(Amplitude,Pattern,:);
            
            PatternX=electrodeMap.getXPosition(Pattern);
            PatternY=electrodeMap.getYPosition(Pattern);
            
            SubplotColumn=(PatternX+945)/120;
            SubplotRow=(PatternY+450)/120;
            
            ConvertedData=ChannelDataToColors2(n,X,Y);
            subplot('position',[SubplotColumn*0.06+0.01 SubplotRow*0.06+0.02 0.054 0.048]);   
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
            %
        end       
    end        
end