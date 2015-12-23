% scripts for plotting moving bar analysis results

clear all
close all

%% parameters

mNo = 37;
ipiAnalysisDone = false;
basePath = '/snle/lab/Experiments/Array/Analysis/2011-06-24-5/';

%% more parameters

%pathToNeuronsFile = '/snle/lab/Experiments/Array/Analysis/2011-06-24-5/data001-from-data000/data001-from-data000.neurons';
pathToEiFile = [basePath 'data006/data006.ei'];

dataNo = 15;
pathToElecResps = [basePath 'data' sprintf('%03.0f', dataNo) '/'];
%cellIDsVisMB =      [2   33   183   393   604   738]; %from moving bar visual stim (data001-from-data000.neurons)
cellIDsVisRF =      [1    33   183   332   214   784]; %from data006.neurons
cellElectrodes =    [60   5    14    21    26    48];
patternNos     =    [1    2    3     4     5     6];
stimAmps =          [0.9  1.9  0.65  0.7   0.8   0.73];
recElecs =          [1    2    13    23    26    50];  %main recording electrodes used in elecResp analysis
isolOccurences =    [82   106  121   112   106   130]; %occurence number (_o[]) corresponding to isolated stimulus
spikeDetThresh =    [0.5  0.5  0.5   0.5   0.5   0.5]; % in fraction of peak signal

nCells = length(cellIDsVisRF);

movieRounding = [1 0; 2 1; 3 2; 4 0; 5 1; 6 2; 7 NaN]; %which rounding interval was used for stimulus starting at this movie (0=no rounding, NaN=isolated pulses)

movToAnalyze = [37 44 51 40 47 54]; %corresponds to 3 highest amplitudes for 1 ms rounding moving bars

%% figure out which o (occurance) numbers correspond with this movie number for each pattern

oNumbers = cell(1,6);

files = dir(pathToElecResps);

%collects all folders corresponding to p1 files
for jj = 1:length(patternNos)
    fileIndex = 0;
    for ii = 1:length(files)
        if strfind(files(ii).name, ['p' num2str(jj) '_o'])==1
            if exist([pathToElecResps filesep files(ii).name filesep files(ii).name '_m' num2str(mNo)],'file')
                fileIndex = fileIndex + 1;
                oNumbers{jj}(fileIndex) = str2double(files(ii).name(strfind(files(ii).name,'o')+1:end));
            end
        end
    end
end
clear fileIndex files

keyboard

%% load elecResps for occurrence corresponding to isolated stim pulses

for ii = 1:length(cellIDsVisRF)
    tmp = load([pathToElecResps filesep 'elecResp_n' num2str(cellIDsVisRF(ii)) '_p' num2str(patternNos(ii))...
        '_o' num2str(isolOccurences(ii)) '.mat']);
    elecRespIso{ii} = tmp.elecResp;
end

%% 


%% load eis

eiFile = edu.ucsc.neurobiology.vision.io.PhysiologicalImagingFile(pathToEiFile);

ei = cell(1,nCells);
for ii = 1:nCells
    ei{ii} = eiFile.getImage(cellIDsVisRF(ii));
    ei{ii} = reshape(ei{ii}(1, 2:end, :), 64, []);
end

%%
%plot overlaid response curves for isolated stim (to choose a good movie
%number)
if 0
    figure; hold on
    cellColors = hsv(length(elecRespIso));
    for ii = 1:length(elecRespIso)
        plot(elecRespIso{ii}.analysis.successRates,...
            'o', 'markerFaceColor', cellColors(ii,:))
    end
end

%% load inter-pulse spike times

%save_path = ['/Analysis/lauren/2011-06-24-5/data' sprintf('%03.0f', dataNo) '-m' num2str(mNo) '-filtered/'];

if ipiAnalysisDone
    load([save_path filesep 'm' num2str(mNo) '_IPIspikes.mat'])
else
    %full_path = ['/snle/lab/Experiments/Array/Analysis/2011-06-24-5/data015-m' num2str(mNo) '-filtered/data' sprintf('%03.0f', dataNo)];
    %save_path = ['/snle/lab/Experiments/Array/Analysis/2011-06-24-5/data015-m' num2str(mNo) '-filtered/'];
    full_path = [basePath 'data' sprintf('%03.0f', dataNo) '-m' num2str(mNo) '-filtered/data' sprintf('%03.0f', dataNo)];
    save_path = [basePath 'data' sprintf('%03.0f', dataNo) '-m' num2str(mNo) '-filtered/'];

    rawFile=edu.ucsc.neurobiology.vision.io.RawDataFile(full_path);
    totalSamples = getRawDataFileSize(full_path);
    
    spikeTimes = cell(1,nCells);
    data = int16(zeros(65, totalSamples));
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
    
    %poor-man's spike finding
    for ii = 1:nCells
        ipiSpikesTmp = ipiSpikeSort(data, ei{ii}(recElecs(ii),:), recElecs(ii), spikeDetThresh(ii));
        keyboard %if you're not happy with the sorting, redo it now!
        
        ipiSpikes(ii).spikeTimes = ipiSpikesTmp.spikeTimes;
        ipiSpikes(ii).cellID = cellIDsVisRF(ii);
        ipiSpikes(ii).recElec = recElecs(ii);
        ipiSpikes(ii).eiPath = pathToEiFile;
        ipiSpikes(ii).spikeDetThresh = spikeDetThresh(ii); %in terms of fraction of peack EI min signal
        ipiSpikes(ii).spikeDetThreshDAQ = ipiSpikesTmp.spikeDetThreshDAQ;
    end
    
    save([save_path filesep 'm' num2str(mNo) '_IPIspikes.mat'], 'ipiSpikes')
end

%% plot rasters
