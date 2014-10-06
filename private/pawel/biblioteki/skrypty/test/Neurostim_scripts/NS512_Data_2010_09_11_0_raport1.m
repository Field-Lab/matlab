NS_GlobalConstants=NS_GenerateGlobalConstants(512);

FolderName='G:\uncompressed\2010-09-14-0\';
DataFilename='008';

full_path='G:\uncompressed\2010-09-14-0\data008';
rawFile=edu.ucsc.neurobiology.vision.io.RawDataFile(full_path);

cd(FolderName);
Movie=82;
%for i=1:224
ChunkData=NS_MovieData(DataFilename,Movie,NS_GlobalConstants);
[PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,MovieData]=NS_DecodeMovieDataChunk(ChunkData);

Channels=[506 507 336 331];
Channels=[457 449 458 465 372 17 364 368];
Data=zeros(length(Channels),20,10000);
for Movie=106:201
    ChunkData=NS_MovieData(DataFilename,Movie,NS_GlobalConstants);
    [PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,MovieData]=NS_DecodeMovieDataChunk(ChunkData);
    for j=1:50
        j
        RawData=int16(rawFile.getData(MovieBegin+RepetPeriod*(j-1),RepetPeriod)');  
        a=RawData(Channels+1,:);
        Data(:,j,:)=a;
    
        FullName=['D:\Home\Pawel\analysis\2010-09-14-0\data008\mv_'  num2str(Movie)];
        %save FullName Data;
        fid=fopen(FullName,'wb','ieee-le');                                    
        fwrite(fid,Data,'int16');
        fclose(fid);        
    end
end
                        
break;

for k=1:length(Channels)
    m=reshape(Data(k,:,:),50,10000);
    plot(t,(m-k*150)');
    hold on;
end