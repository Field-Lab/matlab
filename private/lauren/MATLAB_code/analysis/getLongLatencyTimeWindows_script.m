clear all

%save out .bin file with regions of data between 5 Hz pulses for
%long-latency response analysis (data high-pass filtered to shorten
%artifact regions)

%dataset info

cd /Data/2011-06-24-0/
%cd /tmp/Data/lauren/2011-06-24-0/
FileName = '004';
outpath = ['/snle/lab/Experiments/Array/Analysis/2011-06-24-0/data' FileName '-filtered'];


%parameters of windowing
windowLength = 100*20; %how long of a time window after 5 Hz pulse to extract (including zero-ed out region)
plotLength = 20*50; %length to plot when checking filtering results (doesn't effect saved data)
delayLength = 50; %time after stim pulse to zero out (in samples)


%unchanging parameters
ChipAddresses=[30 31];
NumberOfChannelsPerChip=32;
CurrentRanges=[0.066 0.266 1.07 4.25 16.9 67.1 264 1040];
Fs=20000;
NS_GlobalConstants = struct('SamplingFrequency',Fs,'ChipAddresses',ChipAddresses,'NumberOfChannelsPerChip',NumberOfChannelsPerChip,'CurrentRanges',CurrentRanges);

%extract windows and filter
[filtData movieNumberTimes] = NS_GetLongLatencyTimeWindows(FileName, NS_GlobalConstants, windowLength, delayLength, plotLength);

keyboard

%save filtered data
inpath = [pwd filesep 'data' FileName filesep 'data' FileName '000.bin']; %only use this for header info, so just look at first bin file
%outpath = [pwd filesep 'data' FileName '-filtered'];
outfile = ['data' FileName '000-filtered.bin'];
if ~exist(outpath, 'file')
    mkdir(outpath)
end

save([outpath filesep 'data' FileName '-movieTimes.mat'], 'movieNumberTimes')

fin = fopen(inpath);

header_info = zeros(1,3);
for ii = 1:3
    temp_info = fread(fin, 1, 'uint32', 'b');
    header_info(ii) = temp_info;
end
header_length = header_info(3);


% Rewind the file back to the beginning, read the raw header data
fseek(fin, 0, 'bof');
rawheader = int8(fread(fin, header_length, 'int8'));
fclose(fin);

fout = fopen([outpath filesep outfile], 'w');
fwrite(fout, rawheader, 'int8');
fclose(fout); % Need to do this to get it to flush the output buffer


% append actual data to saved header
write_vision_bin_file(filtData, outpath, outfile, 'append_flag', true)

disp(['based on length of extracted data, there should be ' num2str(size(filtData,2)/(2000*5)) ' repetitions x amplitudes in this dataset'])

% Get data:
%    firstSample = 0;
%    numSamples = 10000;
%    data = rdf.getData(firstSample, numSamples);
%
% Write data back out:
%    outpath = 'Users/peterli/MATLAB'
%    % The data saver uses the header information and decides to save to
%    % /Users/peterli/MATLAB/2011-07-14-6/data001.bin
%    % I didn't mean for it to do this, but seems reasonable.
%
%    edu.ucsc.neurobiology.vision.matlab.Matlab.saveRawData('Users/peterli/MATLAB', rdfh, data);
%
% Check that the data saved out look right:
%    tdf = edu.ucsc.neurobiology.vision.io.RawDataFile('/Users/peterli/MATLAB/2011-07-14-6');
%    all(tdf.getData(0,1) == rdf.getData(0,1))
%    all(tdf.getData(9999,1) == rdf.getData(9999,1))
%    % Both return 1, so the first and 10000th samples match
%
%

        