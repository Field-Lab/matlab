NS_GlobalConstants=NS_GenerateGlobalConstants(61);


cd D:\2008-08-26-0; 
ChunkData=NS_MovieData('008',126,NS_GlobalConstants);
[PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,MovieData]=NS_DecodeMovieDataChunk(ChunkData);
NumberOfEvents=length(MovieData)/3
clear Events;

for k=1:NumberOfEvents
            index=(k-1)*3;
            t=MovieData(index+1);
            PatternNumber=MovieData(index+2);
            Events(k)=PatternNumber;                     
            %Traces(k,j,:,:)=RawData(Channels+1,t+1:t+TraceLength); %always 7 channels
end

plot(Events,'bd-');
grid on
MovieBegin
rawFile=edu.ucsc.neurobiology.vision.io.RawDataFile('D:\2008-08-26-0\data008');
%RawData=int16(rawFile.getData(MovieBegin+RepetPeriod*(j-1),RepetPeriod)');
size(RawData)