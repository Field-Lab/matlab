function [PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,Data]=NS_DecodeMovieDataChunk(ChunkData);

offset=6; %the value 6 come sfrom the fact that 6 values on the beginning on the data array are not for events definition

PDChunkNumber=ChunkData(1); %number of chapter in corresponding SP file
MovieBegin=ChunkData(2)*(2^30)+ChunkData(3); %beginning of first iteration of the movie
RepetNumber=ChunkData(4); %number of repetitions of the movie;
RepetPeriod=ChunkData(5)*(2^30)+ChunkData(6); %period of repetitions in T units (defined in header fo SP file); ignored if defined in movie file header chunk
Data=ChunkData(offset+1:length(ChunkData)); %actual stimulation events, each event having 3 values:
%time in T units, with respect to the current movie chunk
%number of pattern in SP file
%scaling factor