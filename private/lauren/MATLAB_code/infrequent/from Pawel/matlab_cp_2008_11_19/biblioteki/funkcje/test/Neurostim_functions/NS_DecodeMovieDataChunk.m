function [PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,Data]=NS_DecodeMovieDataChunk(ChunkData);

offset=6; %the value 6 come sfrom the fact that 6 values on the beginning on the data array are not for events definition

PDChunkNumber=ChunkData(1);
MovieBegin=ChunkData(2)*(2^30)+ChunkData(3); %beginning of first iteration of the movie
RepetNumber=ChunkData(4); %number of repetitions of the movie;
RepetPeriod=ChunkData(5)*(2^30)+ChunkData(6); 
Data=ChunkData(offset+1:length(ChunkData));