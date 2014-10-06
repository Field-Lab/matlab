cd C:\home\Pawel\nauka\analiza\robocze\2011-07-02-0

Electrodes=[];
for i=1:200
    MovieData=NS_MovieData('001',i,NS_GlobalConstants);
    [PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,MovieData]=NS_DecodeMovieDataChunk(MovieData);
    Electrodes=[Electrodes MovieData(2:3:48)];
end

SE=size(Electrodes);
length(unique(reshape(Electrodes,1,SE(1)*SE(2))))