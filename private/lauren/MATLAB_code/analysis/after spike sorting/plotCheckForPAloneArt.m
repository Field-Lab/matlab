function plotCheckForPAloneArt(pAloneBin, pAmps, neuronID, patternNos, thresholds, threshStds, details, varargin)

% separates primary-alone data into 6 subsets and calculates
% threshold based on subset of data, to see if there is any systematic
% deviation of p-alone from linearity


p = inputParser;

p.addRequired('pAloneBin', @islogical)
p.addRequired('pAmps', @isnumeric)
p.addRequired('neuronID', @isnumeric)
p.addRequired('patternNos', @isnumeric)
p.addRequired('thresholds', @isnumeric)
p.addRequired('threshStds', @isnumeric)
p.addRequired('details', @isstruct)


p.addParamValue('xPlotLim', [-1 1], @isnumeric);
p.addParamValue('yPlotLim', [0 1], @isnumeric);
p.addParamValue('showSecAloneThresh', false, @islogical);
p.addParamValue('fitWithinSecAloneThresh', false, @islogical);


p.parse(pAloneBin, pAmps, neuronID, patternNos, thresholds, threshStds, details, varargin{:})

params = p.Results;


pAloneNeg = find(pAloneBin & pAmps < 0);

%if sum(pAloneBin) ~= 1
if length(pAloneNeg) ~= 1;
    error('More or less than 1 primary-alone pattern detected...error in logic?  SELF-DESTRUCT SEQUENCE....10.....9.....8.....')
end

% extracts data from p-alone elecResp
pAlonePattern = patternNos(pAloneNeg);
temp = load(['elecResp_n' num2str(neuronID) '_p' num2str(pAlonePattern)]);
elecResp = temp.elecResp;

nMovies = length(elecResp.stimInfo.movieNos);

% divides data into 6 subsamples (grouped chronologically)
subSamp = cell(6,1);
subSampLat = cell(6,1);
for jj = 1:6
    subSamp{jj} = [];
end
movieInd = 0;
for ii = 1:nMovies
    if ~isempty(elecResp.analysis.type{ii})
        movieInd = movieInd+1;
        nPerSubSamp = ones(6,1)*floor(elecResp.stimInfo.nPulses(ii)/6);
        %figure out how many pulses are left and distribute them
        remainder = elecResp.stimInfo.nPulses(ii) - nPerSubSamp(1)*6;
        for kk = 1:remainder
            nPerSubSamp(kk) = nPerSubSamp(kk) + 1;
        end
        startInd = 0;
        for jj = 1:6
            subSampLat{jj}{movieInd} = elecResp.analysis.latencies{ii}(startInd+1:startInd+nPerSubSamp(jj));
            subSampTemp(1,1) = abs(elecResp.stimInfo.stimAmps(ii));
            subSampTemp(2,1) = sum(subSampLat{jj}{movieInd}~=0)/nPerSubSamp(jj);
            subSampTemp(3,1) = nPerSubSamp(jj);
            subSamp{jj} = [subSamp{jj} subSampTemp];
            startInd = startInd+nPerSubSamp(jj);
        end
    end
end

subSampThresh = zeros(1,6);
subSampCurveWidths = zeros(1,6);
subSampStDev = zeros(1,6);

% calculate erf fit + bootstrap-based standard deviation for each subsample
for jj = 1:6
    erfParams = erfFitter(subSamp{jj}, 2, -1);
    
    subSampThresh(jj) = -erfParams(2)/erfParams(1);
    subSampCurveWidths(jj) = -sqrt(2)*erfParams(2);
    
    % package data in a usable way for bootstrapThresh (yes this is
    % ridiculous)
    elecRespDummy.stimInfo.movieNos = ones(size(subSamp{jj}, 2),1);
    for ii = 1:size(subSamp{jj}, 2)
        elecRespDummy.analysis.type{ii} = 'blah!';
    end
    elecRespDummy.stimInfo.stimAmps = subSamp{jj}(1,:);
    elecRespDummy.analysis.latencies = subSampLat{jj};
    
    subSampStDev(jj) = bootstrapThresh(elecRespDummy, 100);
end

%plots results of subsampling analysis
threshTmp = thresholds;
threshStdsTmp = threshStds;
detailsTmp = details;

figure('position', [100 100 400 300]*2)
axesPos = {[20 220 160 80]*2, [20 120 160 80]*2, [20 20 160 80]*2, [220 220 160 80]*2, [220 120 160 80]*2, [220 20 160 80]*2};
for jj = 1:6
    axes('units', 'pixels', 'position', axesPos{jj})
    if jj == 1
        title('subsample-based primary-alone thresholds')
    end
    threshTmp(pAloneNeg) = subSampThresh(jj);
    threshStdsTmp(pAloneNeg) = subSampStDev(jj);
    detailsTmp.curveWidths(pAloneNeg) = subSampCurveWidths(jj);
    
    plotLinearityDataVsModel(threshTmp, threshStdsTmp, detailsTmp, jj,...
        'fitWithinSecAloneThresh', params.fitWithinSecAloneThresh, 'showSecAloneThresh', params.showSecAloneThresh)
    set(gca, 'xlim', params.xPlotLim, 'ylim', params.yPlotLim)
end
