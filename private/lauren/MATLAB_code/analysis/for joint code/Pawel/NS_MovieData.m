function MovieData=NS_MovieData(FileName, MovieNumber, varargin)
%This function finds all the time points for which pulse on given electrode
%was generated (started). The function is looking for this time points in
%one movie. 
%Input arguments:
%FileName - in the form of '003';
%-Channel - number of electrode, count from 1 (the function takes care of
%the problem with numbering of channels in raw data files, where actually
%channel number one is the TTL channel; so if you want electrode number
%one, just set this variable to 1);
%-MovieNumber - which movie in given Movie File should be searched;
%-NS_GlobalConstants - the stadard structure.
%Outputs:
%-Timings: just points in time, in reference to the beginning of the run,
%not to the beginning of given movie!
%-PDChunkNumber - in which Pattern Data Chunk in corresponding Pattern File
%the patterns called by given movie are defined.


p = inputParser;

p.addRequired('FileName', @ischar)
p.addRequired('MovieNumber', @isnumeric)

p.addOptional('NS_GlobalConstants', [], @isstruct) %currently unused, but left in for backward compatibility

p.addParamValue('fullPath', [], @ischar) %specifies absolute location of movie file, including file name (FileName is ignored)

p.parse(FileName, MovieNumber, varargin{:})

params = p.Results;

%% 'cleaned' version

if isempty(params.fullPath)
    filename_movie = ['movie' FileName];
else
    filename_movie = params.fullPath;
end
    
fid0=fopen(filename_movie,'r','b');
readMHchunk(fid0); %reads header information


%read in and ignore all the "previous" movies -this is just to find the
%position in the file where the data for THE movie begins. Should be
%changed to "skip" later.
for i=1:MovieNumber-1
    ID=fread(fid0,8,'int8')';
    if ID==[75 116 5 96 -84 122 -59 -64]; %if this is a SC chunk...
        ChunkSize=fread(fid0,1,'int64'); %read in the chunk size
        fread(fid0,ChunkSize,'int32'); %just to move file pointer along
        error('command chunk found in the movie file');
    elseif ID==[114 -69 27 -4 99 66 -12 -123] %MD chunk
        ChunkSize=fread(fid0,1,'int64'); %just to move file pointer along
        fread(fid0,ChunkSize,'int32'); %just to move file pointer along
    end
end

%Now, read THE movie in
ID=fread(fid0,8,'int8')';
if ID==[75 116 5 96 -84 122 -59 -64]; %if this is a SC chunk...
    ChunkSize=fread(fid0,1,'int64'); %just to move file pointer along
    fread(fid0,ChunkSize,'int32'); %just to move file pointer along
    error('command chunk found in the movie file');
elseif ID==[114 -69 27 -4 99 66 -12 -123] %MD chunk
    ChunkSize=fread(fid0,1,'int64');
    MovieData=fread(fid0,ChunkSize,'int32');
end








%% original version
% 
% ChipAddresses=NS_GlobalConstants.ChipAddresses;
% NumberOfChannelsPerChip=NS_GlobalConstants.NumberOfChannelsPerChip;
% CurrentRanges=NS_GlobalConstants.CurrentRanges;
% Fs=NS_GlobalConstants.SamplingFrequency;
% 
% filename_movie=['movie' FileName];
% l=length(filename_movie);
% SPfilename=[filename_movie(1:l-8) 'pattern' filename_movie(l-2:l)];
% 
% fid0=fopen(filename_movie,'r','b');
% %ID=fread(fid0,8,'int8')'
% header=readMHchunk(fid0);
% NumberOfMovies=header.number_of_movies;
% 
% offset=6; %the value 6 come sfrom the fact that 6 values on the beginning on the data array are not for events definition
% 
% %read in and ignore all the "previous" movies -this is just to find the
% %position in the file where the data for THE movie begins. Should be
% %changed to "skip" later.
% for i=1:MovieNumber-1
%     index=0;    
%     ID=fread(fid0,8,'int8')';
%     if ID==[75 116 5 96 -84 122 -59 -64]; %if this is a SC chunk...
%         ChunkSize=fread(fid0,1,'int64'); %read in the chunk size
%         commands=fread(fid0,ChunkSize,'int32');        
%     elseif ID==[114 -69 27 -4 99 66 -12 -123] %MD chunk
%         ChunkSize=fread(fid0,1,'int64');
%         data=fread(fid0,ChunkSize,'int32');
%     end
% end
% 
% %Now, read THE movie in
% ID=fread(fid0,8,'int8')';
% if ID==[75 116 5 96 -84 122 -59 -64]; %if this is a SC chunk...
%     ChunkSize=fread(fid0,1,'int64'); %read in the chunk size
%     commands=fread(fid0,ChunkSize,'int32');        
% elseif ID==[114 -69 27 -4 99 66 -12 -123] %MD chunk
%     ChunkSize=fread(fid0,1,'int64');
%     data=fread(fid0,ChunkSize,'int32');
%     %reading in the movie parameters:
%     PDChunkNumber=data(1);
%     MovieBegin=data(2)*(2^30)+data(3); %beginning of first iteration of the movie
%     RepetNumber=data(4); %number of repetitions of the movie;
%     RepetPeriod=data(5)*(2^30)+data(6); %period of movie repetitions;   
% end
% MovieData=data;