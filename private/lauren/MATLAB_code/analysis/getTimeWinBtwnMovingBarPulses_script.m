clear all

%save out file with containing regions of voltage traces between applied electrical stim pulses

% data is high-pass filtered to shorten artifact regions, and artifact
% regions are zeroes out
%
% data can be saved out to a bin file that can be analyzed by VISION or to
% a .mat file for subsequent MATLAB analysis



delayLength = 30; %time after stim pulse to zero out (in samples); 50 is usually safe
%30 doesn't get rid of artifact on all stim elecs, but is good enough for
%all main rec elecs

% dataset info

cd /tmp/Data/lauren/2011-07-14-0/
FileName = '013';

%cd /tmp/Data/lauren/2013-05-28-3/
%FileName = '013';

%cd /Data/2011-06-24-5/
%FileName = '015';

%parameters of windowing -- unused??? (check NS_getTimeWinBtwnPulses)
plotLength = 20*10; %length to plot when checking filtering results (doesn't affect saved data)
plotElec = [34];

%plotElec = [1 2 14 23 26 50]; %which electrode(s) to plot when checking for remaining artifact--to suppress plotting, set to []

%highest movies are 41, 44
for iMov = [16 23 30 5 12 19 26 33];
    cd /tmp/Data/lauren/2011-07-14-0/

    
    outpath = ['/snle/lab/Experiments/Array/Analysis/2011-07-14-0/data' FileName '-m' num2str(iMov) '-filtered' filesep 'data' FileName];
    outpathMeanSub = [outpath 'meanSub'];
    %outpath = ['/Analysis/2011-06-24-5/data' FileName '-m' num2str(iMov) '-filtered' filesep 'data' FileName];
    %outpath = ['/snle/lab/Experiments/Array/Analysis/2011-06-24-5/data' FileName '-m' num2str(iMov) '-filtered' filesep 'data' FileName];
    
    %extract windows and filter
    [filtData, filtDataMeanSub] = NS_getTimeWinBtwnPulses(FileName, delayLength, iMov, 'plotLength', plotLength, 'plotElec', plotElec);
        
    %save filtered data
    inpath = [pwd filesep 'data' FileName filesep 'data' FileName '000.bin']; %only need this for header info, so just look at first bin file
    %outpath = [pwd filesep 'data' FileName '-filtered'];
    outfile = ['data' FileName '000.bin'];
    if ~exist(outpath, 'file')
        mkdir(outpath)
    end
    if ~exist(outpathMeanSub, 'file')
        mkdir(outpathMeanSub)
    end
    
    %save([outpath filesep 'data' FileName '-movieTimes.mat'], 'movieNumberTimes')
    
    %determine the header length from the original .bin file
    fin = fopen(inpath);
    %header_info = zeros(1,3);
    for ii = 1:3
        temp_info = fread(fin, 1, 'uint32', 'b');
        if ii == 3
            header_length = temp_info;
        end
    end    
    
    % Rewind the file back to the beginning, read the raw header data
    fseek(fin, 0, 'bof');
    rawheader = int8(fread(fin, header_length, 'int8'));
    fclose(fin);
    
    fout = fopen([outpath filesep outfile], 'w');
    fwrite(fout, rawheader, 'int8');
    fclose(fout); % Need to do this to get it to flush the output buffer
    
    s = whos('filtData');
    if s.bytes < 1e8; %under 100 MB
        %append actual data to saved header
        write_vision_bin_file(filtData, outpath, outfile, 'append_flag', true)
    else %might be too big for java heap space (check system limits)
        nPieces = ceil(s.bytes/1e8); %number of pieces required to keep each below 100 MB
        pieceSamples = ceil(size(filtData,2)/nPieces); %number of samples per piece
        for ii = 1:nPieces
            sampStart = (ii-1)*pieceSamples+1;
            sampEnd = ii*pieceSamples;
            if ii < nPieces
                write_vision_bin_file(filtData(:,sampStart:sampEnd), outpath, outfile, 'append_flag', true)
            else
                write_vision_bin_file(filtData(:,sampStart:end), outpath, outfile, 'append_flag', true)
            end
        end
    end
        
    fout = fopen([outpathMeanSub filesep outfile], 'w');
    fwrite(fout, rawheader, 'int8');
    fclose(fout); % Need to do this to get it to flush the output buffer
    
    if s.bytes < 1e8; %under 100 MB
        write_vision_bin_file(filtDataMeanSub, outpathMeanSub, outfile, 'append_flag', true)
    else %might be too big for java heap space (check system limits)
        nPieces = ceil(s.bytes/1e8); %number of pieces required to keep each below 100 MB
        pieceSamples = ceil(size(filtDataMeanSub,2)/nPieces); %number of samples per piece
        for ii = 1:nPieces
            sampStart = (ii-1)*pieceSamples+1;
            sampEnd = ii*pieceSamples;
            if ii < nPieces
                write_vision_bin_file(filtDataMeanSub(:,sampStart:sampEnd), outpathMeanSub, outfile, 'append_flag', true)
            else
                write_vision_bin_file(filtDataMeanSub(:,sampStart:end), outpathMeanSub, outfile, 'append_flag', true)
            end
        end
    end
    
    disp(['extracted ' num2str(size(filtData,2)/(20000)) ' seconds of data'])
    
end


%% sample code provided by Peter (for reference)
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

        