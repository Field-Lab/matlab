function [PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,Data]=NS_DecodeMovieDataChunk(ChunkData)

offset=6; %the value 6 comes from the fact that 6 values on the beginning on the data array are not for events definition

% Points to a location in the patterns file that corresponds to the
% information about this pattern.
PDChunkNumber=ChunkData(1);

MovieBegin=ChunkData(2)*(2^30)+ChunkData(3); %beginning of first iteration of the movie CHUNK
RepetNumber=ChunkData(4); %number of repetitions of the movie chunk
RepetPeriod=ChunkData(5)*(2^30)+ChunkData(6); %period of movie chunk repetitions
Data=ChunkData(offset+1:length(ChunkData)); %the actual movie chunk data (time at which pattern is played, pattern ID, 1)