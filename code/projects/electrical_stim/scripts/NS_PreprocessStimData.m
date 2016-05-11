% This script organizes raw electrical stimulation data according to
% stimulation patterns. A folder for each pattern is created, and data 
% surrounding current pulses for all repetitions are saved in a new file
% for each movie chunk / amplitude
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

javaaddpath('/Volumes/Lab/Development/vision7/Vision.jar');
system = 'stim512';  %stim512 or stim64
%system = 'stim64';
datasets = {'2016-04-21-6'};
datasetsdata = {[2:5 9 10 14 15]};
datasetslen = 1;

for j = 1:datasetslen; 
	% rawDataDir = uigetdir('/Volumes/Data', 'Select raw data directory'); 
	%datNam = '2015-10-29-3';
	datNam = datasets{j};
	rawDataDir = ['/Volumes/Data/' datNam];
	if ~strcmp(rawDataDir(end),filesep)
			rawDataDir = [rawDataDir filesep];
	end

	% Points to the directory of the output.
	%WritePathBase = uigetdir('/Volumes/Analysis', 'Select your output directory'); 
	WritePathBase = ['/Volumes/Analysis/' datNam];
	if ~strcmp(WritePathBase(end),filesep)
			WritePathBase = [WritePathBase filesep];
	end

	% Appends this number to 'data ---'
	%fileNos = 1;
	fileNos = datasetsdata{j};


	% length of trace to save after each pulse (in samples)
	% At 20 samples/millisecond, 100 samples = 5 milliseconds
	traceLength = 100; %changed from 70 on 2010-03-10

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	switch system
			case 'stim512'
					numElectrodes = 512; 
			case 'stim64'
					numElectrodes = 64; 
	end

	% Sets some parameter values depending on what system is in use.
	NS_GlobalConstants = NS_GenerateGlobalConstants(numElectrodes);

	for i = fileNos
			% Set the number of zeros in the file name.
			if i < 10
					FileName = ['00' num2str(i)];
			elseif i < 100
					FileName = ['0' num2str(i)];
			else
					FileName = num2str(i);
			end
			
			rawDataPath = [rawDataDir 'data' FileName];
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
	 
			%%%%%% standard %%%%%%
	%     All repetitions of a particular pattern, at a particular amplitude, inside a particular movie chunk.
	%     For example: Not appropriate for moving bar.
	%     Generates the "p m" files (the pattern-movie files), the
	%     stimulus_files, and the status_files

	NS_SaveRespForMovieAllPatternsAllChannelsNew_sasi(rawDataPath,WritePath,NS_GlobalConstants,traceLength);
			
			%%%%%% special %%%%%%%
			
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
end
