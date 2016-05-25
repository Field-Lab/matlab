NS_GlobalConstants=NS_GenerateGlobalConstants(512);
FileName='I:\analysis\slices\2013-12-12-3-PH\movie001'
Movies=ones(512,450)*1000;
for m=2:451
    m
    ChunkData=NS_MovieData_GlobalPath(FileName,m,NS_GlobalConstants);
    [PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,MovieData]=NS_DecodeMovieDataChunk(ChunkData);
    PatternsUsed=MovieData(2:3:length(MovieData));
    Movies(PatternsUsed,m-1)=m;
end

fid = fopen('D:\Home\Pawel\analysis\slices\2013\2013-12-12-3-PH\data001preproc\MoviesPatterns','wb','ieee-le.l64');
fwrite(fid,min(Movies'),'int32');
fclose(fid);
