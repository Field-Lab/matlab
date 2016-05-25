NS_GlobalConstants=NS_GenerateGlobalConstants(500);
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(500);

NumberOfStimulatedElectrodes=zeros(NumberOfAmplitudes,512);
StimulatedElectrodes=zeros(NumberOfAmplitudes,512,512);
AllPatterns=[]
for Amplitude=25%2:NumberOfAmplitudes
    Amplitude    
    for MovieSeqeunce=1:NumberOfMovieSequences
        Movie=(Amplitude-1)*NumberOfMovieSequences2+MovieSeqeunce
        ChunkData=NS_MovieData_GlobalPath(MovieFile,Movie,NS_GlobalConstants);
        [PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,MovieData]=NS_DecodeMovieDataChunk(ChunkData);  
        
        PatternsForMovie=unique(MovieData(2:3:length(MovieData)));
        AllPatterns=[AllPatterns' PatternsForMovie']';
        for p=1:length(PatternsForMovie)            
            Pattern=PatternsForMovie(p)                        
            [DataTraces,ArtifactDataTraces,Channels]=NS_ReadPreprocessedData(DataPath,ArtifactDataPath,0,Pattern,Movie,1500,0);            
            
            Events=NS512_DetectSpikes(DataTraces,threshold,threshold/2);            
            [a,b]=find(Events>0);            
            n=hist(b,[1:512]);            
            
            %ElectrodesWithSpikes=find(n>25);
            %NumberOfStimulatedElectrodes(Amplitude,Pattern)=length(ElectrodesWithSpikes);
            StimulatedElectrodes(Amplitude,Pattern,:)=n;                    
        end               
    end        
end