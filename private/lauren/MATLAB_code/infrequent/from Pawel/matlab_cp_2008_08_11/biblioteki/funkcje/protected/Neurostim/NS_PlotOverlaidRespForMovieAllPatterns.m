function [PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,y]=NS_PlotOverlaidRespForMovieAllPatterns(FileName,WritePath,FigureProperties,NS_GlobalConstants);
% This function generates quickly plots of overlaid responses for all
% patterns and all movies FOR 1-EL. STIMULATION!!!! It is assumed that the
% number of pattern called in the movie is the number of stimulating
% electrode.
ArrayID=1;

ChipAddresses=NS_GlobalConstants.ChipAddresses;
NumberOfChannelsPerChip=NS_GlobalConstants.NumberOfChannelsPerChip;
CurrentRanges=NS_GlobalConstants.CurrentRanges;
Fs=NS_GlobalConstants.SamplingFrequency;

full_path=[pwd '/' 'data' FileName];
rawFile=edu.ucsc.neurobiology.vision.io.RawDataFile(full_path);
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(1);

TraceLength=60;

filename_movie=['movie' FileName];
l=length(filename_movie);
SPfilename=[filename_movie(1:l-8) 'pattern' filename_movie(l-2:l)];

fid0=fopen(filename_movie,'r','b'); %opens stimulation movie file
%ID=fread(fid0,8,'int8')'

% returns a struct containing 'number_of_movies' (total number of Movie Data (MD) chunks in the file)
% 'period' (period of repetition for each movie, set to 0 if periods are to be defined
% individually, otherwise the period value defined in the movie data chunk is ignored)
% (defined in T units, which are defined in header of corresponding Pattern File)
header=readMHchunk(fid0);
NumberOfMovies=header.number_of_movies;

for i=1:NumberOfMovies %unless specific movie numbers defined!!
    ID=fread(fid0,8,'int8')'; %reads the first 8 int8 values
    if ID==[75 116 5 96 -84 122 -59 -64]; %if this is a SC (stimchip commands) chunk...
        ChunkSize=fread(fid0,1,'int64'); %read in the chunk size (in bytes)
        commands=fread(fid0,ChunkSize,'int32'); %reads in the rest of the SC chunk
    elseif ID==[114 -69 27 -4 99 66 -12 -123] %MD (movie data) chunk
        ChunkSize=fread(fid0,1,'int64'); %read in the chunk size (in bytes)
        ChunkData=fread(fid0,ChunkSize,'int32'); %reads in the rest of the MD chunk
        
        %reading in the movie parameters:
        %PDChunkNumber = number of chapter in corresponding SP file
        %MovieBegin = beginning of first iteration of the movie in T units
        %RepetNumber = number of repetitions in T units
        %MovieData = actual stimulation events, each event having 3 values:
        %time in T units, with respect to the current movie chunk
        %number of pattern in SP file
        %scaling factor
        [PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,MovieData]=NS_DecodeMovieDataChunk(ChunkData);
    end
    MovieBegin;
    NumberOfEvents=length(MovieData)/3;
    
    [Patterns,PatternsIndexes,Status]=ReadPatternDataChunk(SPfilename,PDChunkNumber,NS_GlobalConstants);
    NumberOfPatterns=length(PatternsIndexes)
    ChannelsToShow=ones(NumberOfPatterns,7)*9;
    Traces=zeros(NumberOfPatterns,RepetNumber,7,TraceLength);  % always 7 electrodes will be shown - should be more flexible !!!!!!
                 
    for j=1:RepetNumber
        RawData=double(rawFile.getData(MovieBegin+RepetPeriod*(j-1),RepetPeriod)');
        %size(RawData)
        for k=1:NumberOfEvents
            index=(k-1)*3;
            t=MovieData(index+1);
            PatternNumber=MovieData(index+2);
            
            [Name,Channels]=NS_PatternAmplitudes(Patterns,PatternsIndexes,Status,PatternNumber,NS_GlobalConstants);
            Names{PatternNumber}=Name;
                        
            neighbors=electrodeMap.getAdjacentsTo(PatternNumber,1); % * * * * GOOD ONLY FOR 1-EL. STIMULUS!! might be less than 7
            ChannelsToShow(PatternNumber,1:length(neighbors))=neighbors';   
            Channels=ChannelsToShow(PatternNumber,:);  % always 7 channels, because such is the size of "ChannelsToShow" array
            Traces(PatternNumber,j,:,:)=RawData(Channels+1,t+1:t+TraceLength); %always 7 channels
        end
    end
    
    for l=1:NumberOfPatterns % this loop should be running only over patterns used in this movie !!!!
        Channels=electrodeMap.getAdjacentsTo(l,1);  % * * * * GOOD ONLY FOR 1-EL. STIMULUS!! might be less than 7
        TracesToShow=reshape(Traces(l,:,1:length(Channels),:),RepetNumber,length(Channels),TraceLength);% might be <7         
        y=NS_PlotManySignaturesOnArrayLayoutNew(TracesToShow,Channels,ones(1,RepetNumber),ArrayID,FigureProperties,NS_GlobalConstants);        
        hj=gcf;
        set(hj, 'PaperOrientation', 'portrait');
        
        FullName=[WritePath '\' 'p' num2str(l) '_m' num2str(i) Names{l}]
        print(hj, '-dtiff', FullName);
    end    
end  