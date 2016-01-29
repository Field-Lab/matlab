function [PDChunkNumber,movieBegin,numRepetitions,repetitionPeriod,Data] = getMovieData(movieFileName, movieChunkNumber, ~)
%GETMOVIEDATA finds all the time points in one specified movie at which a 
% a pulse was generated on a given electrode. 

% This function combines NS_MovieData and NS_DecodeMovieData LGrosberg

%Input arguments:
% movieFileName: full path to the labview output movie file, 'movieXXX';
% movieChunkNumber: specifies the movie file chunk (i.e., the sequence of patterns applied during the specified refresh period) in a given movie file 'movieXXX' to decode
% NS_GlobalConstants: Optional parameter, unused  but defined for
% compatibility with previous versions of this function. 

%Outputs:
% PDChunkNumber: in which Pattern Data Chunk in corresponding Pattern File
% the patterns called by given movie are defined.
% movieBegin: beginning of the first iteration of the movie CHUNK
% repetitionPeriod: period of movie chunk repetitions
% Data:  is structured as: [ x    y    z ...]'
%               x = event time (in samples)
%                   y = ID of the pattern applied
%                       z = the number 1 (not used)
% Every 3 entries define the information for 1 event.

    
fid0=fopen(movieFileName,'r','b');
readMHchunk(fid0); 
 
% Move file pointer to the correct "movie file chunk" location in the movie file 
for i=1:movieChunkNumber-1
    ID=fread(fid0,8,'int8')';
    if ID==[75 116 5 96 -84 122 -59 -64]; %if this is a SC chunk...
        ChunkSize=fread(fid0,1,'int64'); %read in the chunk size
        fread(fid0,ChunkSize,'int32'); %just to move file pointer along
        error('command chunk found in the movie file');
    elseif ID==[114 -69 27 -4 99 66 -12 -123] % MD chunk
        ChunkSize=fread(fid0,1,'int64'); %just to move file pointer along
        fread(fid0,ChunkSize,'int32'); %just to move file pointer along
    end
end

% Read the current movie file chunk
ID=fread(fid0,8,'int8')';
if ID==[75 116 5 96 -84 122 -59 -64]; %if this is a SC chunk...
    ChunkSize=fread(fid0,1,'int64'); %just to move file pointer along
    fread(fid0,ChunkSize,'int32'); %just to move file pointer along
    error('command chunk found in the movie file');
elseif ID==[114 -69 27 -4 99 66 -12 -123] %MD chunk
    ChunkSize=fread(fid0,1,'int64');
    ChunkData=fread(fid0,ChunkSize,'int32');
end


offset = 6; % the first 6 values of the movie chunk data array do not define events
PDChunkNumber = ChunkData(1);
movieBegin = ChunkData(2)*(2^30)+ChunkData(3); 
numRepetitions = ChunkData(4); 
repetitionPeriod = ChunkData(5)*(2^30)+ChunkData(6); 
Data = ChunkData(offset+1:length(ChunkData)); %the actual movie chunk data (time at which pattern is played, pattern ID, 1) 

end

