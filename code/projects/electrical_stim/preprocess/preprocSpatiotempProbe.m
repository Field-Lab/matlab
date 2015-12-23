function [PDChunkNumber, MovieBegin, nRepeats, repeatPeriod]=preprocSpatiotempProbe(FileName, WritePath, NS_GlobalConstants,...
    traceLength, centerChannels, stimPeriod)

% preprocesses raw electrical stimulation (STIM64) data for spatiotemporal probe
%
% trace segments begin with pre-pulse
% pattern number in file names are replaced with p(attern)[number of stim elec]_pre[number of
% prepulse elec]_d[delay between pulse onsets, in samples]
%
% warning: this code relies heavily on expected timing of pulses (stim
% pulse always occurs on multiple of stimPeriod/nCenters)
%
% base function was originally written by Pawel, but has been modified by
% Lauren in 2010-03

%ArrayID=1;

currentRanges = NS_GlobalConstants.CurrentRanges;

% The location of the data (set to the current directory).
full_path=[pwd filesep 'data' FileName] %#ok<NOPRT>

% Creates a java object (rawFile) that let's you interact with the raw data files (sort of the stitched together bin files).
rawFile=edu.ucsc.neurobiology.vision.io.RawDataFile(full_path);

% Returns the number of samples contained in the raw data file.
totalSamples = getRawDataFileSize(full_path);

% Generates the full filenames for movie and pattern.
filename_movie=['movie' FileName] %#ok<NOPRT>
SPfilename=['pattern' FileName] %#ok<NOPRT>

% Open the movie file.
fid0=fopen(filename_movie,'r','b');
% Read through the header.
header=readMHchunk(fid0);
% Number of movie chunks x number of amplitudes.
nMovies=header.number_of_movies;

% Calculates the number of channels, 1:64.
Channels=1:64;
% Initialize a stopwatch.
t0 = clock;

% Converts stimPeriod (the length between stimulation trials) from ms to samples. (20 samples/ms)
stimPeriod = stimPeriod*20;
% Calculates the number of center stimulation electrodes.
nCenters = length(centerChannels);

% Iterates through all the movies.
% It always saves an extra movie at the end that should be discarded.

% Unless specific movie numbers defined!!
% Minus 1 because the last movie is expected to be empty.
for i=1:nMovies-1 
    
    % Displays progress/time left.
    if i>1
        finished = (i-1)/(nMovies-1); % proportion of files created so far
        disp(sprintf('finished writing %0.1f%% of files', finished*100))
        tnow = clock;
        timeElapsed = etime(tnow, t0); % time elapsed since loop started
        estimatedTimeLeft = (timeElapsed/finished)*(1-finished);
        disp(sprintf('estimated time left: %0.1f seconds',estimatedTimeLeft))
    end
        
    % Read the chunk type ID.
    ID=fread(fid0,8,'int8')';
    
    % Checks for MD chunks by returning an error if ID is an SC chunk.
    if ID==[75 116 5 96 -84 122 -59 -64]; %#ok<BDSCA,BDSCI>
        error('command chunk found in the movie file');
        %ChunkSize = fread(fid0,1,'int64'); %read in the chunk size
        %commands=fread(fid0,ChunkSize,'int32');        
    elseif ID==[114 -69 27 -4 99 66 -12 -123] %#ok<BDSCA,BDSCI> %MD chunk
        
        % need to leave these in so that the file pointer is in the correct position for the next
        % iteration
        ChunkSize=fread(fid0, 1, 'int64'); %size of MD chunk
        fread(fid0, ChunkSize, 'int32'); %just to move file pointer along
        
        %reading in the movie parameters:
        ChunkData=NS_MovieData(FileName,i,NS_GlobalConstants); %stimulus information: times at which patterns are played
        
        % interprets first 6 values in chunkdata
        % see NS_DecodeMovieDataChunk for descriptions
        [PDChunkNumber,MovieBegin,nRepeats,repeatPeriod,MovieData]=NS_DecodeMovieDataChunk(ChunkData);

    end
    nEvents=length(MovieData)/3; %number of pattern applications in movie
    
    [patterns,patternsIndexes,Status]=ReadPatternDataChunk(SPfilename,PDChunkNumber,NS_GlobalConstants); %#ok<NASGU>
    FullName_status=[WritePath filesep 'status_files' filesep 'status_m' num2str(i)];
    save(FullName_status,'Status')
    
    % Calculate the length of traces that will be saved for each pattern
    % (or set of patterns).
    TraceLengthToSave = MovieData(4) - MovieData(1);
    if TraceLengthToSave >= stimPeriod
        
        TraceLengthToSave = traceLength;
                
    else
        delayInd = 0;
        stimDelays = zeros(1,nEvents/2);
    
        for ii = 1:6:length(MovieData)
            delayInd = delayInd+1;
            stimDelays(delayInd) = MovieData(ii+3) - MovieData(ii);
        end

        if any(stimDelays ~= mean(stimDelays)) 
            error('stimDelays are non-uniform in the movie');
        end

        
        TraceLengthToSave = TraceLengthToSave + traceLength;
    end
    

    
    Traces=zeros(nEvents, nRepeats, length(Channels), TraceLengthToSave, 'int16'); 
    
    for j=1:nRepeats %iterates through number of times movie is played
        
        %checks to make sure raw data file has all of the samples that the
        %stimulus files thinks it does...
        if totalSamples <= MovieBegin+repeatPeriod*j+TraceLengthToSave
            warndlg(['Raw data file ' full_path ' does not have as many samples as expected by stimulus files']);
            return
        end
        
        % added traceLength samples to retrieved portion of data to prevent index
        % out-of-bounds error
        RawData=int16(rawFile.getData(MovieBegin + repeatPeriod*(j-1), repeatPeriod + TraceLengthToSave)');
        subgroupInfo = struct([]);
        iSubgroup = 0;
        
        for k = 1:nEvents
                
            index=(k-1)*3;
            t = MovieData(index+1);
            patternNumber = MovieData(index+2);
            
            
            %get timing of preceding pulse
            if k > 1
                tPrev = MovieData(index-2);
            else
                tPrev = [];
            end
            
            %get timing of following pulse
            if k < nEvents
                tNext = MovieData(index+4);
            else
                tNext = [];
            end
            
            % Pattern1 will contain information about the stimulus for this
            % pattern.
            pattern1=NS_ExtractPatternFromPatternDataChunk(patterns,patternsIndexes,patternNumber);
            % determines which channels were used in stimulus
            stimChannelsP1 = [];
            for m = 1:65
                if ~any(m==[1, 10, 26, 58]) && any(pattern1(m).data(1,:)) %legitimate channel with stim
                    stimChannelsP1 = [stimChannelsP1 pattern1(m).channel]; %#ok<AGROW>
                end
            end
            
            if isempty(stimChannelsP1) %empty pattern
                error('encountered an empty pattern - check that preprocessing code will still work before proceeding!')
                
            elseif length(stimChannelsP1) == 2 %simultaneous stim
                %check to make sure timing of pulse pair is a multiple of
                %stimPeriod/nCenters as expected
                if mod(t, stimPeriod/nCenters) ~= 0
                    error('simultaneous stimulation occurs at an unexpected time')
                end
                
                if nCenters == 1
                    centerChannelBin = (stimChannelsP1 == centerChannels(1));
                else
                    centerChannelBin = (stimChannelsP1 == centerChannels(1)) | (stimChannelsP1 == centerChannels(2));
                end
                centerChannel = stimChannelsP1(centerChannelBin);
                otherChannel = stimChannelsP1(~centerChannelBin);
                
                
                currentRange = currentRanges(Status.ChannelsStatus(otherChannel).range + 1)/127;
                preAmps = pattern1(otherChannel+1).data(1,:).*pattern1(otherChannel+1).data(3,:)*currentRange;
                preAmp = max(abs(preAmps));
                preAmpStr = num2str(preAmp, '%0.3f');
                
                preAmpStr(strfind(preAmpStr, '.')) = '_';
                
                patternName = ['p' num2str(centerChannel) '_pre' num2str(otherChannel) '_d0_a' num2str(preAmpStr)];
                patternLength = traceLength;
                
                %check for already-named subgroup or create new one
                foundMatch = 0;
                for m = 1:length(subgroupInfo)
                    if strcmpi(patternName, subgroupInfo(m).patternName)
                        subgroupInfo(m).eventBin(k) = 1;
                        foundMatch = 1;
                    end
                end
                if ~foundMatch
                    iSubgroup = iSubgroup + 1;
                    subgroupInfo(iSubgroup).eventBin = false(nEvents, 1);
                    subgroupInfo(iSubgroup).eventBin(k) = 1;
                    subgroupInfo(iSubgroup).patternName = patternName;
                    subgroupInfo(iSubgroup).patternLength = patternLength;
                    subgroupInfo(iSubgroup).pattern = pattern1;
                end
                
                skip = 0;

            % prepulse because it doesn't occur at a multiple of stimPeriod/nCenters
            elseif length(stimChannelsP1)==1 && mod(t, stimPeriod/nCenters) ~= 0
                if k == nEvents % current pulse is last pulse in movie
                    error('there is no stimulus pulse following apparent prepulse')
                
                elseif nCenters == 1 && tNext - t > stimPeriod*0.95 % time between current and next pulse is more than 95% of the period
                    error(['there is no stimulus pulse following apparent prepulse within ' num2str(stimPeriod*0.95/20) ' ms, which is 0.95x the full period'])
                elseif nCenters == 2 && tNext - t > stimPeriod*0.45 % time between current and next pulse is more than 45% of the period
                    error(['there is no stimulus pulse following apparent prepulse within ' num2str(stimPeriod*0.45/20) ' ms, which is 0.45x the full period (with 2 center elecs)'])
                elseif mod(tNext, stimPeriod/nCenters) ~= 0
                    error('stimulus pulse following apparent prepulse doesn''t occur at expected time')
                end
                
                %extract pattern data for following (stim) pulse
                pattern2=NS_ExtractPatternFromPatternDataChunk(patterns,patternsIndexes,MovieData(index+5)); %pattern for stim pulse
                stimChannelP2 = [];
                for m = 1:65
                    if ~any(m==[1, 10, 26, 58]) && any(pattern2(m).data(1,:)) %legitimate channel with stim
                        stimChannelP2 = [stimChannelP2 pattern2(m).channel]; %#ok<AGROW>
                    end
                    delay = tNext-t; %between prepulse and stim pulse, in samples
                end
                
                %check to make sure following (stim) pulse involves a
                %single, center electrode
                if ~any(centerChannels == stimChannelP2) || length(stimChannelP2) ~= 1
                    error('unexpected sequence of pulses--aborting')
                end
                

                %determine amplitude of prepulse for pattern naming purposes
                currentRange = currentRanges(Status.ChannelsStatus(stimChannelsP1).range + 1)/127;
                preAmps = pattern1(stimChannelsP1+1).data(1,:).*pattern1(stimChannelsP1+1).data(3,:)*currentRange;
                preAmp = max(abs(preAmps));
                preAmpStr = num2str(preAmp, '%0.3f');
                preAmpStr(strfind(preAmpStr, '.')) = '_';

                
                patternName = ['p' num2str(stimChannelP2) '_pre' num2str(stimChannelsP1) '_d' num2str(delay) '_a' num2str(preAmpStr)];
                
                patternLength = delay + traceLength;
                
                %check for already-named subgroup or create new one
                foundMatch = 0;
                for m = 1:length(subgroupInfo)
                    if strcmpi(patternName, subgroupInfo(m).patternName)
                        subgroupInfo(m).eventBin(k) = 1;
                        foundMatch = 1;
                    end
                end
                if ~foundMatch
                    iSubgroup = iSubgroup + 1;
                    subgroupInfo(iSubgroup).eventBin = false(nEvents, 1);
                    subgroupInfo(iSubgroup).eventBin(k) = 1;
                    subgroupInfo(iSubgroup).patternName = patternName;
                    subgroupInfo(iSubgroup).patternLength = patternLength;
                    
                    %combine pattern data of prepulse and stimpulse
                    pattern = struct;
                    maxPatternLength = 0;
                    for m = 1:65
                        maxPatternLength = max(maxPatternLength, size(pattern2(m).data, 2));
                    end
                    for m = 1:65
                        pattern(m).data = zeros(5, delay+maxPatternLength);%padded by 1 sample before and 2 samples after
                        pattern(m).data(:, 1:size(pattern1(m).data,2)) = pattern1(m).data;
                        pattern(m).data(:, delay+1:delay+size(pattern2(m).data,2)) = pattern2(m).data;
                        pattern(m).channel = pattern1(m).channel;
                    end
                    subgroupInfo(iSubgroup).pattern = pattern;
                end
                
                skip = 0;
                           
            %center alone (w/o prepulse) because current and previous pulses occur at a multiple of stimPeriod/nCenters, and current
            %stimulus includes single, center electrode
            elseif length(stimChannelsP1) == 1 && any(centerChannels == stimChannelsP1) && mod(t, stimPeriod/nCenters) == 0 && (k==1 || mod(tPrev, stimPeriod/nCenters) == 0) %center alone
                if k > 1
                    prevDelay = t-tPrev;
                    if nCenters == 1 && prevDelay < stimPeriod*0.95 %less than 0.95*full period since last pulse
                        error(['center alone stim expected to be more than ' num2str(stimPeriod*0.95/20) ' ms from previous pulse'])
                    elseif nCenters == 2 && prevDelay < stimPeriod*0.45 %less than 0.45*full period since last pulse
                        error(['center alone stim expected to be more than ' num2str(stimPeriod*0.45/20) ' ms from previous pulse'])
                    end
                end
                
                patternName = ['p' num2str(stimChannelsP1)];
                patternLength = traceLength;
                
                %check for already-named subgroup or create new one
                foundMatch = 0;
                for m = 1:length(subgroupInfo)
                    if strcmpi(patternName, subgroupInfo(m).patternName)
                        subgroupInfo(m).eventBin(k) = 1;
                        foundMatch = 1;
                    end
                end
                if ~foundMatch
                    iSubgroup = iSubgroup + 1;
                    subgroupInfo(iSubgroup).eventBin = false(nEvents, 1);
                    subgroupInfo(iSubgroup).eventBin(k) = 1;
                    subgroupInfo(iSubgroup).patternName = patternName;
                    subgroupInfo(iSubgroup).patternLength = patternLength;
                    subgroupInfo(iSubgroup).pattern = pattern1;
                end
                
                skip = 0;
                
            %stim pulse following pre-pulse
            elseif (length(stimChannelsP1) == 1 && any(centerChannels == stimChannelsP1)) && mod(t, stimPeriod/nCenters) == 0 
                if k == 1
                    error('no prepulse prior to apparent stim pulse')
                elseif nCenters == 1 && t-tPrev > stimPeriod*0.95 %more than 0.95*full period since last pulse
                    error(['no prepulse within ' num2str(stimPeriod*0.95/20) ' ms prior to apparent stim pulse'])
                elseif nCenters == 2 && t-tPrev > stimPeriod*0.45 %more than 0.45*full period since last pulse
                    error(['no prepulse within ' num2str(stimPeriod*0.45/20) ' ms prior to apparent stim pulse'])
                end
                skip = 1;
            else
                keyboard
                error('unexpected sequence of pulses')
            end
            
            if ~skip
                Traces(k,j,:,1:patternLength)=RawData(Channels+1, t+1:t+patternLength);
            end
        end
    end
  
    
    
    for j = 1:length(subgroupInfo)
        pattern=subgroupInfo(j).pattern; %#ok<NASGU>

        TracesToSave=Traces(subgroupInfo(j).eventBin,:,1:length(Channels),1:subgroupInfo(j).patternLength);
        TracesToSave = reshape(TracesToSave, nRepeats*sum(subgroupInfo(j).eventBin),length(Channels),subgroupInfo(j).patternLength);
        
        STTS=size(TracesToSave);

        a=reshape(TracesToSave,STTS(1)*STTS(2)*STTS(3),1);
        b=zeros(1000,1);
        b(1)=STTS(1);
        b(2)=STTS(2);
        b(3)=STTS(3);
        b(3+1:3+length(Channels))=Channels';
        o=[b' a'];

        if ~exist([WritePath filesep subgroupInfo(j).patternName], 'file')
            mkdir([WritePath filesep subgroupInfo(j).patternName])
        end
        
        FullName=[WritePath filesep subgroupInfo(j).patternName filesep subgroupInfo(j).patternName '_m' num2str(i)];
        fid=fopen(FullName,'wb','ieee-le');
        fwrite(fid,o,'int16');
        fclose(fid);

        FullName_pattern=[WritePath filesep 'pattern_files' filesep 'pattern' subgroupInfo(j).patternName(2:end) '_m' num2str(i)];
        save(FullName_pattern,'pattern');
    end 
    clear Traces;
    clear TracesToSave;
end  