function [PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,ClusterFileName]=NS_PreprocessDataNew512v6main519(pathToRawData,WritePath,WritePathFigs,ArrayID,FigureProperties,NS_GlobalConstants,makePlots,PatternsToRead,ChannelsToRead,MovieNumberOffset);
% This function generates quickly plots of overlaid responses for all
% patterns and all movies FOR 1-EL. STIMULATION!!!! It is assumed that the
% number of pattern called in the movie is the number of stimulating
% electrode.
% PatternsToRead - if movie includes 64 patterns, but we want to process
% only some of them, the indexes should be specifiec here. The aray of
% indexes is not array of patterns really, but indexes for the array ofPatternsUsed
% patterns that will be loaded by this function from the stimulus files.


rawFile=edu.ucsc.neurobiology.vision.io.RawDataFile(pathToRawData);
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(ArrayID);

TraceLength=140;

[parentstr,datarun,~] = fileparts(pathToRawData); 
datarunNumber = datarun(5:end); 
altFileOrg = fullfile(parentstr,'Electrical', 'Output');
if ~isempty(dir(fullfile(parentstr,'movie*')))
    movieFileName= fullfile(parentstr,['movie' datarunNumber]);
    SPfilename    = fullfile(parentstr,['pattern' datarunNumber]);
elseif ~isempty(dir(fullfile(altFileOrg,'movie*')))
    movieFileName= fullfile(altFileOrg, ['movie' datarunNumber]);
    SPfilename    = fullfile(altFileOrg,['pattern' datarunNumber]);
else
    % User selects directory with the labview output files
    correctFileOrg = uigetdir(parentstr, 'Select directory containing LabVIEW output files ''movie*'' and ''pattern*'''); 
    if correctFileOrg
        movieFileName = fullfile(correctFileOrg,['movie' datarunNumber]);
        SPfilename     = fullfile(correctFileOrg,['pattern' datarunNumber]);
    else
        err = MException('MATLAB:missingLabVIEWFiles', ...
            'Cannot find labview output files that define the times and patterns of electrical stimulation');
        throw(err);
    end
end


fid0=fopen(movieFileName,'r','b');
header=readMHchunk(fid0);
NumberOfMovies=header.number_of_movies

Channels=ChannelsToRead;
t0 = clock;
NumberOfPatternsMax=0;
NumberOfRepetitionsMax=0;

for ii= 1:1:NumberOfMovies-1 %unless specific movie numbers defined!! minus 1 because the last movie is expected to be empty. This loop repeats once for each amplitude (in case of scan data).    
    %a ) estimate how much time is left to complete the loop
  
    if ii>1 
        finished = (ii-1)/(NumberOfMovies-1); % proportion of files created so far
        fprintf('finished writing %0.1f%% of files\n', finished*100);
        tnow = clock;
        timeElapsed = etime(tnow, t0); %time elapsed since loop started
        estimatedTimeLeft = (timeElapsed/finished)*(1-finished);
        fprintf('estimated time left: %0.1f minutes\n',estimatedTimeLeft/60)
    end

    %b) read in single movie_data_chunk    
    ID=fread(fid0,8,'int8')';
    if ID==[75 116 5 96 -84 122 -59 -64]; %if this is a SC (command) chunk...
        error('command chunk found in the movie file');
        ChunkSize=fread(fid0,1,'int64'); %read in the chunk size
        commands=fread(fid0,ChunkSize,'int32');        
    elseif ID==[114 -69 27 -4 99 66 -12 -123] %MD (movie data) chunk
        ChunkSize=fread(fid0,1,'int64');
        ChunkData=fread(fid0,ChunkSize,'int32');
        %reading in the movie parameters:
        ChunkData=NS_MovieData(datarunNumber,ii,NS_GlobalConstants);
        [PDChunkNumber,MovieBegin,RepetNumber0,RepetPeriod,MovieData]=NS_DecodeMovieDataChunk(ChunkData);        
        RepetNumber=min(RepetNumber0,25);
    end    

    %calculate number of events related to specified patterns
    PatternsUsed=MovieData(2:3:length(MovieData));
    
    PatternsIndexesToRead=[]; % 
    for i1=PatternsToRead
        PatternIndexes=find(PatternsUsed==i1);  %number of events corresponding to given pattern
        PatternsIndexesToRead=[PatternsIndexesToRead' PatternIndexes']'; 
    end
              
    NumberOfEvents=length(PatternsIndexesToRead);  % this will not work when some patterns are used more than once in the movie!!!
    Events=zeros(NumberOfEvents,1);
    
    %c) read in corresponding pattern_data_chunk, save status into file, predefine array for RAW data  
  
    [Patterns,PatternsIndexes,Status]=ReadPatternDataChunk(SPfilename,PDChunkNumber,NS_GlobalConstants); %#ok<ASGLU>
    disp(['ReadPatternDataChunk done for Chunk' num2str(PDChunkNumber) ' ' SPfilename]);
    FullName_status=fullfile(WritePath,'status_files', ['status_m' num2str(ii+MovieNumberOffset)]);
    if ~exist(fileparts(FullName_status),'dir')
        mkdir(fileparts(FullName_status))
    end
    save(FullName_status,'Status');
            
    NumberOfPatterns=length(PatternsIndexes);
    NumberOfPatternsMax=max(NumberOfPatternsMax,NumberOfPatterns);
    length(ChannelsToRead);
    numberofbytes = NumberOfEvents*RepetNumber*length(ChannelsToRead)*TraceLength*2 ;   
    Traces=[];
    if numberofbytes>0
        if numberofbytes>100000000 % poniewa? maksymalna wielkosc tablicy to 943 MB, wi?c dla duzych tablic trzeab stworzyc najpierw ma?a tablic?, skonwertowa? na typ int16 i dopiero z takiej sk?ada? du?? talbic? int16
            TracesSmall=int16(zeros(NumberOfEvents,1,length(ChannelsToRead),TraceLength)); 
            Traces=TracesSmall;
            for ble=2:RepetNumber               
                Traces=[Traces TracesSmall];                
            end
        else
            Traces=int16(zeros(NumberOfEvents,RepetNumber,length(ChannelsToRead),TraceLength));  
        end                                
    size(Traces);
    %d) read in RAW data for each repetition, for each stimulation event
    %identify pattern number, and save RAW data coresponding to each event into
    %Traces array. We assume each repetition of the movie includes identical
    %collection of patterns (logical, isn't it).

    for j=1:RepetNumber % one iteration of this loop takes less than 0.5s for 64-channel data (122 events)                
        %RepIndex=j
        RawData=int16(rawFile.getData(MovieBegin+RepetPeriod*(j-1),RepetPeriod+TraceLength)');                        
        for k=1:NumberOfEvents % later on, single iteration should read and organize data for all ... of given pattern in this movie at once (PH, 2010-08-26)
            PatternNumber=PatternsUsed(PatternsIndexesToRead(k));
            t=MovieData(1+(PatternsIndexesToRead(k)-1)*3)  ;                                          
            Events(k)=PatternNumber;
            size(RawData);
            size(RawData(Channels+1,t+1:t+TraceLength));            
            Traces(k,j,:,:)=RawData(Channels+1,t+1:t+TraceLength); %this takes less than 0.1ms for 64-channel data                                 
        end        
    end    
    
    DifferentPatterns=unique(Events);
        
    clear RawData;
    
    for l=1:length(DifferentPatterns) %1:NumberOfPatterns % this loop should be running only over patterns used in this movie !!!!
        ThisPattern=DifferentPatterns(l);  
        WhichEvents=find(Events==ThisPattern); %which events in this movie corresponded to given pattern
        NumberOfRepetitionsMax=max(NumberOfRepetitionsMax,RepetNumber*length(WhichEvents));
        if length(WhichEvents)>0            
            Pattern=NS_ExtractPatternFromPatternDataChunk(Patterns,PatternsIndexes,ThisPattern);                                    
                                    
            TracesToSave=reshape(Traces(WhichEvents,:,1:length(Channels),:),RepetNumber*length(WhichEvents),length(Channels),TraceLength);              
            STTS=size(TracesToSave);  
            if ArrayID==1502
                a=importdata(fullfile(matlab_code_path,'code','projects','electrical_stim','stimuli','519_30um_stimulation','519mess.txt'));
                TracesToSave519=zeros(STTS(1),519,STTS(3));
                TracesToSave519(:,a,:)=TracesToSave;
                clear TracesToSave;
                TracesToSave=TracesToSave519;
                clear TracesToSave519;
                STTS=size(TracesToSave);
            end            
        
            a=reshape(TracesToSave,STTS(1)*STTS(2)*STTS(3),1);
            
            b=zeros(1000,1);
            b(1)=STTS(1);
            b(2)=STTS(2);
            b(3)=STTS(3);
            b(3+1:3+length(Channels))=Channels';
            b(4+length(Channels))=length(WhichEvents);
            o=[b' a'];            
                        
            if ArrayID==1502
                FullName=fullfile(WritePath,['p' num2str(NS512_519Inverse(ThisPattern))],['p' num2str(NS512_519Inverse(ThisPattern)) '_m' num2str(ii+MovieNumberOffset)]);
            else
                FullName=[WritePath filesep 'p' num2str(ThisPattern) filesep 'p' num2str(ThisPattern) '_m' num2str(ii+MovieNumberOffset)];
            end
            if ~exist(fileparts(FullName),'dir')
                mkdir(fileparts(FullName));
            end
            fid=fopen(FullName,'wb','ieee-le');                                 
            fwrite(fid,o,'int16');
            fclose(fid);                    
                        
            if ArrayID==1502
                FullName_pattern= fullfile(WritePath, 'pattern_files',['pattern' num2str(NS512_519Inverse(ThisPattern)) '_m' num2str(ii+MovieNumberOffset)]);
            else
                FullName_pattern=[WritePath filesep 'pattern' num2str(ThisPattern) '_m' num2str(ii+MovieNumberOffset)];  
            end   
            if ~exist(fileparts(FullName_pattern),'dir')
                mkdir(fileparts(FullName_pattern))
            end
            save(FullName_pattern,'Pattern');
            
            % part below - only for 1-el. scan!!
            if makePlots==1
                ChannelsToPlot=electrodeMap.getAdjacentsTo(ThisPattern,1);                
                TracesToShow=TracesToSave(:,ChannelsToPlot,:);
                %c=int16(mean(TracesToShow));
                %for ci=1:STTS(1)
                %    TracesToShow(ci,:,:)=TracesToShow(ci,:,:)-c;
                %end
            
                y = NS_PlotManySignaturesOnArrayLayoutNew(TracesToShow,ChannelsToPlot,ones(1,STTS(1)),ArrayID,FigureProperties,NS_GlobalConstants);
                hj=gcf;              
                set(hj,'PaperUnits','inches');
                set(hj,'PaperSize',[13 9]);
                set(hj,'PaperPosition',[0 0 13 9]);
                FullName=[WritePathFigs '\' 'p' num2str(ThisPattern) '_m' num2str(ii)];
                print(hj, '-dtiff', '-r120', FullName);
            end
        end
    end                                                                          
    clear Traces;
    clear TracesToSave;
    end
end  
fclose(fid0);
% Create cluster file:
ClusterFileName=NS_CreateClusterFile(WritePath,datarunNumber,NumberOfMovies,NumberOfPatternsMax,NumberOfRepetitionsMax);