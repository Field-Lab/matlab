% script to load the filtered data into Matlab and then run it through
% ipiSpikeSort to generate *_IPIspikes.mat file of spike times
%
% assumes standard naming/directory system for filtered data
%   e.g. .../2011-06-24-5/data015-m1-filtered/data015/data015000.bin
%
% saves output file into same directory containing filtered data directory as
% 'm[number]_IPIspikes.mat'


clear all

useMeanSub = false;


%where filtered data exists and where ipiSpikes will be saved
BasePath = '/snle/lab/Experiments/Array/Analysis/2011-07-14-0/';
EstimData = 'data013';

movieChunk = 30; % 30,33,37,40

pathToEiFile = [BasePath 'data005/data005.ei']; %ei file used in elecResp files

% ordering of the following must be consistent 
% (and potentially must be in the order defined by the corresponding 
% pattern numbers, but this SHOULD be checked when combining ipi spikes with elecResp spikes)
cellIDs =           [48   187   304   407   498   618   783   918]; %should match cells in specified ei file
recElecs =          [4    14    24    33    37    45    53    61];  %main recording electrodes used in elecResp analysis
spikeDetThresh =    [0.5  0.5   0.5   0.5   0.5   0.5   0.5   0.5]; % in fraction of peak signal




% %where filtered data exists and where ipiSpikes will be saved
% BasePath = '/snle/lab/Experiments/Array/Analysis/2013-05-28-3/';
% EstimData = 'data013'; %data012 or data013
% movieChunk = 41; %41 or 44
% 
% pathToEiFile = [BasePath 'data010/data010.ei']; %ei file used in elecResp files
% 
% % ordering of the following must be consistent 
% % (and potentially must be in the order defined by the corresponding 
% % pattern numbers, but this SHOULD be checked when combining ipi spikes with elecResp spikes)
% cellIDs =           [166   138   349   437   784   916]; %should match cells in specified ei file
% recElecs =          [6     13    22    31    47    62];  %main recording electrodes used in elecResp analysis
% spikeDetThresh =    [0.5   0.5   0.5   0.5   0.5   0.5]; % in fraction of peak signal




% %where filtered data exists and where ipiSpikes will be saved
% BasePath = '/snle/lab/Experiments/Array/Analysis/2011-06-24-5/';
% EstimData = 'data015';
% 
% movieChunk = 37;
% 
% pathToEiFile = [BasePath 'data006/data006.ei']; %ei file used in elecResp files
% 
% % ordering of the following must be consistent 
% % (and potentially must be in the order defined by the corresponding 
% % pattern numbers, but this SHOULD be checked when combining ipi spikes with elecResp spikes)
% cellIDs =           [1    33   183   332 data010  214   784]; %should match cells in specified ei file
% recElecs =          [1    2    13    23    26    50];  %main recording electrodes used in elecResp analysis
% spikeDetThresh =    [0.5  0.5  0.5   0.5   0.5   0.5]; % in fraction of peak signal

%%

nCells = length(cellIDs);

% load up ei file (used as templates for setting detection thresholds and
% comparing expected spike waveforms to spikes identified in the data)
eiFile = edu.ucsc.neurobiology.vision.io.PhysiologicalImagingFile(pathToEiFile);
ei = cell(1,nCells);
for ii = 1:nCells
    ei{ii} = eiFile.getImage(cellIDs(1,ii));
    ei{ii} = reshape(ei{ii}(1, 2:end, :), 64, []);
end
eiFile.close();

% load the filtered data into Matlab's memory
% have to do it a little bit at a time because of java heap space limitations
if useMeanSub
    pathToFilteredData = [BasePath EstimData '-m' num2str(movieChunk) '-filtered/' EstimData 'meanSub' filesep];
else
    pathToFilteredData = [BasePath EstimData '-m' num2str(movieChunk) '-filtered/' EstimData filesep];
end

rawFile=edu.ucsc.neurobiology.vision.io.RawDataFile(pathToFilteredData);
totalSamples = getRawDataFileSize(pathToFilteredData);

data = int16(zeros(65, totalSamples));
disp(['loading ' num2str(round(totalSamples/20000)) 's of data from ' pathToFilteredData])
for ii = 1:ceil(totalSamples/20000) %load in 1 second at a time to avoid running out of java heap space
    fprintf('%c', '.')
    if mod(ii,100)==0
        fprintf('\n')
    end
    iSamp = (ii-1)*20000;

    if ii <  ceil(totalSamples/20000)
        data(:,iSamp+1:iSamp+20000) = int16(rawFile.getData(iSamp, 20000))';
    else
        data(:,iSamp+1:end) = int16(rawFile.getData(iSamp, totalSamples - iSamp))';
    end
end

%throw out the trigger channel
data(1,:) = [];

% run filtered data through ipiSpikeSort to extract spike times
for ii = 1:nCells
    ipiSpikesTmp = ipiSpikeSort(data, ei{ii}, recElecs(ii), spikeDetThresh(ii));
    
    disp('here''s your chance to redo that cell if you didn''t sort the spikes correctly!')
    keyboard
    
    ipiSpikes(ii).spikeTimes = ipiSpikesTmp.spikeTimes;
    ipiSpikes(ii).cellID = cellIDs(ii);
    ipiSpikes(ii).recElec = recElecs(ii);
    ipiSpikes(ii).eiPath = pathToEiFile;
    ipiSpikes(ii).ei = ei{ii};
    ipiSpikes(ii).spikeDetThresh = spikeDetThresh(ii);
    ipiSpikes(ii).spikeDetThreshDAQ = ipiSpikesTmp.spikeDetThreshDAQ;
end

% save to disk
if useMeanSub
    saveName = [BasePath EstimData '-m' num2str(movieChunk) '-filtered/m' num2str(movieChunk) 'meanSub_IPIspikes.mat'];
else
    saveName = [BasePath EstimData '-m' num2str(movieChunk) '-filtered/m' num2str(movieChunk) '_IPIspikes.mat'];
end

if exist(saveName, 'file')
    button = questdlg(['A file named ' saveName ' already exists!' 10 'Overwrite?'], '','yes', 'no', 'no');
    if strcmp(button,'yes')
        save(saveName, 'ipiSpikes')
    end
else
    save(saveName, 'ipiSpikes')
end