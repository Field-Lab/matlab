function [PDChunkNumber, MovieBegin, nRepeats, repeatPeriod]=NS_SaveRespSpatiotemporalStim(FileName, WritePath, NS_GlobalConstants, traceLength, stimPeriod)

% preprocesses raw electrical stimulation (STIM64) data
%
% this function was originally written by Pawel, but has been modified (in mostly cosmetic ways) by
% Lauren in 2009-09

%
%  This function assumes that all stimuli except for 2nd-final pulses in a
%  sequence occur at times = multiples of stimPeriod
%


full_path=[pwd filesep 'data' FileName] %#ok<NOPRT>

rawFile=edu.ucsc.neurobiology.vision.io.RawDataFile(full_path);


totalSamples = getRawDataFileSize(full_path);

filename_movie=['movie' FileName] %#ok<NOPRT>
l=length(filename_movie);
SPfilename=[filename_movie(1:l-8) 'pattern' filename_movie(l-2:l)] %#ok<NOPRT>

fid0=fopen(filename_movie,'r','b');
header=readMHchunk(fid0);
NumberOfMovies=header.number_of_movies;

Channels=1:64;
t0 = clock;

stimPeriod = stimPeriod*20; %converts from ms to samples

sequenceOrders = {}; %cell array storing all sequence orders used in stimulus

% iterates through stimulus amplitudes ('movies')
for i=1:NumberOfMovies-1 %unless specific movie numbers defined!! minus 1 because the last movie is expected to be empty
    if i>1
        finished = (i-1)/(NumberOfMovies-1); % proportion of files created so far
        disp(sprintf('finished writing %0.1f%% of files', finished*100))
        tnow = clock;
        timeElapsed = etime(tnow, t0); %time elapsed since loop started
        estimatedTimeLeft = (timeElapsed/finished)*(1-finished);
        disp(sprintf('estimated time left: %0.1f seconds',estimatedTimeLeft))
    end
        
    ID=fread(fid0,8,'int8')';
    
    if ID==[75 116 5 96 -84 122 -59 -64]; %if this is a SC chunk...
        error('command chunk found in the movie file');

    elseif ID==[114 -69 27 -4 99 66 -12 -123] %MD chunk
        
        % need to leave these in so that the file pointer is in the correct position for the next iteration
        ChunkSize=fread(fid0, 1, 'int64'); %size of MD chunk
        fread(fid0, ChunkSize, 'int32'); %just to move file pointer along
        
        %reading in the movie parameters:
        ChunkData=NS_MovieData(FileName,i,NS_GlobalConstants); %stimulus information: times at which patterns are played
        
        
        % interprets first 6 values in chunkdata
        % see NS_DecodeMovieDataChunk for descriptions
        [PDChunkNumber,MovieBegin,nRepeats,repeatPeriod,MovieData]=NS_DecodeMovieDataChunk(ChunkData);
    end
    nEvents=length(MovieData)/3; %number of pattern applications in movie
    
    [Patterns,PatternsIndexes,Status]=ReadPatternDataChunk(SPfilename,PDChunkNumber,NS_GlobalConstants); %#ok<NASGU>
    FullName_status=[WritePath filesep 'status_files' filesep 'status_m' num2str(i)];
    save(FullName_status,'Status')
        
    Events = MovieData(2:3:end); %pattern numbers
    EventTimes = MovieData(1:3:end);
    
    
    % determines which stimuli are part of a sequence, which are
    % simultaneous, and which are isolated pulses
    stimTypes = cell(1, nEvents);
    orderNumbers = zeros(1, nEvents); %if part of a sequence, stimulus is assigned a number corresponding to which order it is
    patternChannels = zeros(1, nEvents); %for single-elec stimulus (sequence or isolated), store what channel is used in stimulus
    orderDelays = zeros(1, nEvents);
    sequenceStarts = false(1, nEvents); %true if this pattern is the first in a sequence
    for k = 1:nEvents
        
        %index=(k-1)*3;
        patternNumber = Events(k);
        t = EventTimes(k);
        
        %indeces of patterns that are applied at a multiple of the
        %stimPeriod
        periodMult = find(mod(EventTimes, ones(size(EventTimes))*stimPeriod) == 0);
        
        %determine number of electrodes involved in stimulus
        channelsTemp = [];
        patternTemp=NS_ExtractPatternFromPatternDataChunk(Patterns,PatternsIndexes,patternNumber);
        for mm = 1:65
            if ~any(mm==[1, 10, 26, 58]) && any(patternTemp(mm).data(1,:)) %legitimate channel with stim
                channelsTemp = [channelsTemp patternTemp(mm).channel];
            end
        end
        nChannels = length(channelsTemp);
        clear patternTemp;
        

        if nChannels > 1
            stimTypes{k} = 'simultaneous';
        elseif k ~= nEvents
            %if current or following stimuli don't occur on a multiple of
            %stimPeriod, must be part of a sequence
            if ~any(periodMult == k) || ~any(periodMult == k+1)
                stimTypes{k} = 'sequence';
            else
                stimTypes{k} = 'isolated';
            end
        else %k == nEvents
            if ~any(periodMult == k)
                stimTypes{k} = 'sequence';
            else
                stimTypes{k} = 'isolated';
            end
        end
        
        if strcmpi(stimTypes{k}, 'isolated') || strcmpi(stimTypes{k}, 'sequence')
           patternChannels(k) = channelsTemp;
        end
        clear channelsTemp
            
        
        if strcmpi(stimTypes{k}, 'sequence')
            %figure out the order of electrodes in the sequence
            sequenceStart = periodMult(find(periodMult<=k, 1, 'last'));
            sequenceEnd = periodMult(find(periodMult>k, 1, 'first'))-1;
            iInSequence = sequenceStart:sequenceEnd;
            nInSequence = length(iInSequence);
            
            if sequenceStart == k
                sequenceStarts(k) = true;
            end
            
            sequenceOrderTemp = [];
            for ii = 1:length(iInSequence)
                %determine electrode involved in each pattern
                patternTemp=NS_ExtractPatternFromPatternDataChunk(Patterns,PatternsIndexes,Events(iInSequence(ii)));
                channelsTemp = [];
                for mm = 1:65
                    if ~any(mm==[1, 10, 26, 58]) && any(patternTemp(mm).data(1,:)) %legitimate channel with stim
                        channelsTemp = [channelsTemp patternTemp(mm).channel]; %#ok<AGROW>
                    end
                end
                clear patternTemp
                
                if length(channelsTemp) > 1
                    error('unexpected sequence of stimuli')
                end
                
                sequenceOrderTemp(ii) = channelsTemp; %#ok<AGROW>
            end
            
            %check if this sequence has already been found
            for ii = 1:length(sequenceOrders)
                if all(sequenceOrders{ii} == sequenceOrderTemp)
                    orderNumbers(k) = ii;
                    break
                end
            end
            if orderNumbers(k) == 0 %no match found
                orderNumbers(k) = length(sequenceOrders)+1;
                sequenceOrders{orderNumbers(k)} = sequenceOrderTemp; %#ok<AGROW>
            end
            orderDelays(k) = EventTimes(iInSequence(2)) - EventTimes(iInSequence(1));
        end
    end
    
    fullSeqTraceLength = traceLength+max(orderDelays)*nInSequence;
            
    Traces=int16(zeros(nEvents, nRepeats, 64, traceLength));
    TracesFullSeq = int16(zeros(sum(sequenceStarts), nRepeats, 64, fullSeqTraceLength)); %special traces that contain entire sequences
    
    for j=1:nRepeats %iterates through number of times movie is played
        iSequence = 0;
        %checks to make sure raw data file has all of the samples that the
        %stimulus files thinks it does...
        if totalSamples <= MovieBegin+repeatPeriod*j+traceLength
            warndlg(['Raw data file ' full_path ' does not have as many samples as expected by stimulus files']);
            return
        end
        

        % %% hack!!! (for spatial stim patterns that had incorrect movie lengths specified in stim files)
        %         %hack to fix incorrect MovieBegin and repeatPeriod
        %         MovieBegin = 286000+(i-1)*1300000;
        %
        %         repeatPeriod = 26000;
        %
        % %% end hack
       
        %retrieve windows of data
        RawData=int16(rawFile.getData(MovieBegin + repeatPeriod*(j-1), repeatPeriod + traceLength)');
        for k = 1:nEvents
            t = EventTimes(k);
            Traces(k,j,:,:)=RawData(Channels+1, t+1:t+traceLength);
            
            if sequenceStarts(k)
                iSequence = iSequence+1;
                TracesFullSeq(iSequence, j, :, :) = RawData(Channels+1, t+1:t+fullSeqTraceLength);
            end
        end
    end
            
    
    iSequence = 0;
    for jj=1:nEvents
        Pattern=NS_ExtractPatternFromPatternDataChunk(Patterns,PatternsIndexes,Events(jj));%#ok<NASGU>

        TracesToSave=reshape(Traces(jj,:,1:length(Channels),:),nRepeats,length(Channels),traceLength);
        STTS=size(TracesToSave);

        a=reshape(TracesToSave,STTS(1)*STTS(2)*STTS(3),1);
        b=zeros(1000,1);
        b(1)=STTS(1);
        b(2)=STTS(2);
        b(3)=STTS(3);
        b(3+1:3+length(Channels))=Channels';
        o=[b' a'];


        if strcmpi(stimTypes{jj}, 'simultaneous')
            if ~exist([WritePath filesep 'psim'], 'file')
                mkdir([WritePath filesep 'psim'])
            end
            FullName = [WritePath filesep 'psim' filesep 'psim_m' num2str(i)];
            FullName_pattern=[WritePath filesep 'pattern_files' filesep 'patternsim' '_m' num2str(i)];
        elseif strcmpi(stimTypes{jj}, 'isolated')
            if ~exist([WritePath filesep 'p' num2str(patternChannels(jj))], 'file')
                mkdir([WritePath filesep 'p' num2str(patternChannels(jj))])
            end
            FullName = [WritePath filesep 'p' num2str(patternChannels(jj)) filesep...
                'p' num2str(patternChannels(jj)) '_m' num2str(i)];
            FullName_pattern=[WritePath filesep 'pattern_files' filesep...
                'pattern' num2str(patternChannels(jj)) '_m' num2str(i)];
        elseif strcmpi(stimTypes{jj}, 'sequence')
            if ~exist([WritePath filesep 'p' num2str(patternChannels(jj)) '_d' num2str(orderDelays(jj)) '_o' num2str(orderNumbers(jj))], 'file')
                mkdir([WritePath filesep 'p' num2str(patternChannels(jj)) '_d' num2str(orderDelays(jj)) '_o' num2str(orderNumbers(jj))])
            end
            FullName = [WritePath filesep 'p' num2str(patternChannels(jj)) '_d' num2str(orderDelays(jj)) '_o' num2str(orderNumbers(jj)) filesep...
                'p' num2str(patternChannels(jj)) '_d' num2str(orderDelays(jj)) '_o' num2str(orderNumbers(jj)) '_m' num2str(i)];
            FullName_pattern=[WritePath filesep 'pattern_files' filesep...
                'pattern' num2str(patternChannels(jj)) '_d' num2str(orderDelays(jj)) '_o' num2str(orderNumbers(jj)) '_m' num2str(i)];
        else
            error('unexpected stimType--error in preprocessing code')
        end

       
        fid=fopen(FullName,'wb','ieee-le');
        fwrite(fid,o,'int16');
        fclose(fid);
        clear TracesToSave;
        
        save(FullName_pattern,'Pattern');
        
        %separate set of preproc files to store data for full seqeunces
        if sequenceStarts(jj)
            iSequence = iSequence+1;
            thisSeqTraceLength = traceLength+orderDelays(jj)*nInSequence;
            TracesToSave=reshape(TracesFullSeq(iSequence,:,1:length(Channels),1:thisSeqTraceLength),nRepeats,length(Channels),thisSeqTraceLength);
            STTS=size(TracesToSave);
            
            a=reshape(TracesToSave,STTS(1)*STTS(2)*STTS(3),1);
            b=zeros(1000,1);
            b(1)=STTS(1);
            b(2)=STTS(2);
            b(3)=STTS(3);
            b(3+1:3+length(Channels))=Channels';
            o=[b' a'];
            
            if ~exist([WritePath filesep 'pfullseq_d' num2str(orderDelays(jj)) '_o' num2str(orderNumbers(jj))], 'file')
                mkdir([WritePath filesep 'pfullseq_d' num2str(orderDelays(jj)) '_o' num2str(orderNumbers(jj))])
                
            end
            FullName = [WritePath filesep 'pfullseq_d' num2str(orderDelays(jj)) '_o' num2str(orderNumbers(jj)) filesep...
                'pfullseq_d' num2str(orderDelays(jj)) '_o' num2str(orderNumbers(jj)) '_m' num2str(i)];
            fid=fopen(FullName,'wb','ieee-le');
            fwrite(fid,o,'int16');
            fclose(fid);
        end
        
    end
    clear Traces TracesFullSeq TracesToSave;
end  