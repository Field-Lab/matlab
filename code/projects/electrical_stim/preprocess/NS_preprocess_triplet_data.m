% Extracts all the data recorded in a pulse triplet experiment 

% Extracts the data surrounding current pulses and saves to a new file.
% Does all repetitions of one specific stimulus pulse in a given movie chunk at once.

system = 'stim512'; %'stim64'; %stim512 or stim64
%system = 'stim64';

% Points to the directory of the raw data.
% cd /Volumes/Acquisition/Data/2014-04-15-4
% cd /Volumes/Data/2014-04-10-0
% cd /Volumes/Data/2014-04-15-4
% cd /Volumes/Data/2014-06-04-3
% cd /Volumes/Data/2014-07-08-3
% cd /Volumes/Data/2012-09-24-3/
% cd /Volumes/Acquisition/Data/2014-07-21-test
% cd /Volumes/Data/2014-07-24-0
% cd /Volumes/Data/2012-09-24-0/
cd /Volumes/Data/2015-03-09-0

% Points to the directory of the output.
%WritePathBase = '/Analysis/2012-09-18-2/';
%WritePathBase = '/braid/snle/data/2011-06-24-5/';
% WritePathBase = '/Volumes/Analysis/2014-04-15-4/';
% WritePathBase = '/Volumes/Analysis/2014-06-04-3/';
% WritePathBase = '/Volumes/Analysis/2014-07-08-3/';
% WritePathBase = '/Volumes/Analysis/delete/2012-09-24-3/'; 
% WritePathBase = '/Volumes/Acquisition/Analysis/2014-07-21-test/';
% WritePathBase = '/Volumes/Analysis/2014-07-24-0/long-chunks/';
% WritePathBase = '/Volumes/Analysis/2012-09-24-0/';
WritePathBase = '/Volumes/Analysis/2015-03-09-0/'; 

% Appends this number to 'data ---'
fileNos = [14];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sets some parameter values depending on what system is in use.
if strcmpi(system, 'stim512')
    % ChipAddresses is an arbitrary designation number for each chip.
    ChipAddresses=31:-1:24;
    NumberOfChannelsPerChip=64;
    % Each array has a particular ArrayID.
    ArrayID=500;
elseif strcmpi(system, 'stim64')
    ChipAddresses=[30 31];
    NumberOfChannelsPerChip=32;
    ArrayID=1;
else
    error('no other systems currently supported')
end

% These 3 parameters never change - correspond to hardware.
% Upper limit of current ranges. Constant within a movie chunk (for an
% electrode?)
CurrentRanges=[0.066 0.266 1.07 4.25 16.9 67.1 264 1040];
% Fs = sampling rate (in Hz)
Fs=20000;

NS_GlobalConstants=struct('SamplingFrequency',Fs,'ChipAddresses',ChipAddresses,...
    'NumberOfChannelsPerChip',NumberOfChannelsPerChip,'CurrentRanges',CurrentRanges);


for i = fileNos
    % Set the number of zeros in the file name.
    if i < 10
        FileName = ['00' num2str(i)];
    elseif i < 100
        FileName = ['0' num2str(i)];
    else
        FileName = num2str(i);
    end

    WritePath = [WritePathBase 'data' FileName];
    
    % Determine if WritePath exists, if not then make it.
    if ~exist(WritePath, 'file')
        mkdir(WritePath)
        disp(['Created new directory: ' WritePath])
    end

    if ~exist([WritePath filesep 'pattern_files'], 'file')
        mkdir([WritePath filesep 'pattern_files'])
        disp(['Created new directory: ' WritePath filesep 'pattern_files'])
    end

    if ~exist([WritePath filesep 'status_files'], 'file')
        mkdir([WritePath filesep 'status_files'])
        disp(['Created new directory: ' WritePath filesep 'status_files'])
    end

    % length of trace to save after each pulse (in samples)
    % At 20 samples/millisecond, 100 samples = 5 milliseconds
    traceLength = 600; %30 ms of data

    
    %%%%%% standard %%%%%%
%     All repetitions of a particular pattern, at a particular amplitude, inside a particular movie chunk.
%     For example: Not appropriate for moving bar.
%     Generates the "p m" files (the pattern-movie files), the
%     stimulus_files, and the status_files

% NS_SaveRespForMovieAllPatternsAllChannelsNew(FileName,WritePath,NS_GlobalConstants,traceLength);
    
    %%%%%% special %%%%%%%
[PDChunkNumber, MovieBegin, nRepeats, repeatPeriod]=NS_SaveRespForMovieTriplets(FileName, WritePath,...
    NS_GlobalConstants, traceLength);
    % Used if it is inappropriate to group together all the patterns in a
    % given movie chunk, for example, if the same pattern is used in one
    % movie chunk multiple times in a moving bar.  Generates the
    % distinction of 'occurrences' within a repetition.
%NS_SaveRespEachOccurrenceSeparately(FileName,WritePath,NS_GlobalConstants, traceLength);
    
    % Used if a stimulus is interrupted in the middle of running by
    % shifting the baseline movie number by 'mShift'
    %mShift = 88;
%NS_SaveRespForMovieAllPatternsAllChannelsNew_shifted_m(FileName, WritePath, NS_GlobalConstants, traceLength, mShift)

    % Detects the frequency of stimulation and outputs as part of the file
    % name.
%NS_SaveRespFrequencyScan(FileName,WritePath,NS_GlobalConstants, traceLength);

    % Was used for spatiotemporal sequences, 'stimPeriod' is time between
    % each sequence applied.
    %stimPeriod = 50; %in ms
%NS_SaveRespSpatiotemporalStim(FileName, WritePath, NS_GlobalConstants, traceLength, stimPeriod);
   
    % Used for spatiotemporal probe on 'centerChannel', with 'stimPeriod'
    % between each pulse/pre-pulse pair.
%     centerChannels = 22;
%     stimPeriod = 128; %in ms
% preprocSpatiotempProbe(FileName, WritePath, NS_GlobalConstants, traceLength, centerChannels, stimPeriod);
%     
    % Use for extracting the timing of occurrences from moving bar data.
    % Patterns:     [1  2   3   4   5   6]
    % Occurrences:  [13 14  21  18  17  21]  
%     totalPatterns = 6;
%     numberOfOccurrences = 21;
%     [PDChunkNumber, MovieBegin, nRepeats, repeatPeriod, pattern_application_times] =...
%         extract_stim_occurrence_timing(FileName, WritePath, NS_GlobalConstants, totalPatterns);
%     
    % Something.
    %fakePreprocessStimData(FileName, WritePath, 100, [1 10000])
end