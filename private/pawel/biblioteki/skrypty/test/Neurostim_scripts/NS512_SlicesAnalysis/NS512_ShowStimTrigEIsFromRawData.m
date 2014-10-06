full_path='D:\Home\Data\slices\2010-09-14-0\data002'; %define path to raw data file
rawFile=edu.ucsc.neurobiology.vision.io.RawDataFile(full_path); 

NS_GlobalConstants=NS_GenerateGlobalConstants(500);
MovieFilePath='D:\Home\Data\slices\2010-09-14-0\movie002';
patternID=30;

movieIDList = NS512_GetMoviesWithPatternIDNew(MovieFilePath,patternID);
WhichMovie=10;

MovieData4=NS_MovieData_GlobalPath(MovieFilePath,movieIDList(WhichMovie),NS_GlobalConstants);
[PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,MovieData]=NS_DecodeMovieDataChunk(MovieData4);
SMD=size(MovieData);
PatternNumbers=MovieData(2:3:SMD(1));
WhichPatterns=find(PatternNumbers==patternID);

WaveformLength=600;
spikes3=zeros(RepetNumber,512,WaveformLength);
for i=1:length(WhichPatterns)
    TimeIndex=MovieData(1+(WhichPatterns(i)-1)*3);
    PulseTimes=MovieBegin+TimeIndex+(0:RepetNumber-1)*RepetPeriod;
    
    for t=1:RepetNumber
        PulseTime=PulseTimes(t);
        d0=rawFile.getData(PulseTime,WaveformLength)';        
        spikes3(t,:,:)=d0(2:513,:);        
    end    
end