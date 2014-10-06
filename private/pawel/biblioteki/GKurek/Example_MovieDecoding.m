cd F:\2010-08-02-1
MovieData4=NS_MovieData('003',4,NS_GlobalConstants);
[PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,MovieData]=NS_DecodeMovieDataChunk(MovieData4);
PDChunkNumber
MovieBegin
RepetNumber
RepetPeriod

fid0=fopen('movie003','r','b');
%ID=fread(fid0,8,'int8')'
header=readMHchunk(fid0);
NumberOfMovies=header.number_of_movies
fclose(fid0);