NS_GlobalConstants=NS_GenerateGlobalConstants(512);

Channels=[20 124 122 135 185 188 217 246 292 290 308 318 330 372 384 428 433 444 433 445 446 449 445 455 1 457 459 457 456 467 453 497];
Ch1=20;
Ch2=446;
ChannelsPlot1=electrodeMap.getAdjacentsTo(Ch1,1)';
ChannelsPlot2=electrodeMap.getAdjacentsTo(Ch2,1)';
Channels=[ChannelsPlot1 ChannelsPlot2];
Movies=[82:8:136];
ilosc=50;

FolderName='G:\uncompressed\2010-09-14-0\';

DataFilename='002';
full_path='G:\uncompressed\2010-09-14-0\data002';
rawFile=edu.ucsc.neurobiology.vision.io.RawDataFile(full_path);
cd(FolderName);
Data=zeros(length(Channels),20,10000);
for Movie=Movies%1:136
    ChunkData=NS_MovieData(DataFilename,Movie,NS_GlobalConstants);
    [PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,MovieData]=NS_DecodeMovieDataChunk(ChunkData);
    for j=1:ilosc
        j
        RawData=int16(rawFile.getData(MovieBegin+RepetPeriod*(j-1),RepetPeriod)');  
        a=RawData(Channels+1,:);
        Data(:,j,:)=a;
    
        FullName=['D:\Home\Pawel\analysis\2010-09-14-0\data002new\mv_'  num2str(Movie)];
        fid=fopen(FullName,'wb','ieee-le');                                    
        fwrite(fid,Data,'int16');
        fclose(fid);        
    end
end


DataFilename='005';
full_path='G:\uncompressed\2010-09-14-0\data005';
rawFile=edu.ucsc.neurobiology.vision.io.RawDataFile(full_path);
cd(FolderName);
Data=zeros(length(Channels),25,10000);
for Movie=Movies%1:136
    ChunkData=NS_MovieData(DataFilename,Movie,NS_GlobalConstants);
    [PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,MovieData]=NS_DecodeMovieDataChunk(ChunkData);
    for j=1:ilosc
        j
        RawData=int16(rawFile.getData(MovieBegin+RepetPeriod*(j-1),RepetPeriod)');  
        a=RawData(Channels+1,:);
        Data(:,j,:)=a;
    
        FullName=['D:\Home\Pawel\analysis\2010-09-14-0\data005new\mv_'  num2str(Movie)];
        fid=fopen(FullName,'wb','ieee-le');                                    
        fwrite(fid,Data,'int16');
        fclose(fid);        
    end
end


DataFilename='008';
full_path='G:\uncompressed\2010-09-14-0\data008';
rawFile=edu.ucsc.neurobiology.vision.io.RawDataFile(full_path);
cd(FolderName);
Data=zeros(length(Channels),20,10000);
for Movie=Movies%1:200
    ChunkData=NS_MovieData(DataFilename,Movie,NS_GlobalConstants);
    [PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,MovieData]=NS_DecodeMovieDataChunk(ChunkData);
    for j=1:ilosc
        j
        RawData=int16(rawFile.getData(MovieBegin+RepetPeriod*(j-1),RepetPeriod)');  
        a=RawData(Channels+1,:);
        Data(:,j,:)=a;
    
        FullName=['D:\Home\Pawel\analysis\2010-09-14-0\data008new\mv_'  num2str(Movie)];
        fid=fopen(FullName,'wb','ieee-le');                                    
        fwrite(fid,Data,'int16');
        fclose(fid);        
    end
end