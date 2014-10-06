function [latencies threshold thresholdStd] = fitErfToSubset(elecRespPath, subsetBin)

% !!!!!!!!!!!! THIS IS AN OLD FUNCTION THAT IS ONLY USED BY OUTDATED CODE (spatioTempProbeFigures_old) !!!!!!!!!!!!!!!!!!!!!!!
% script to fit curve to subset of trials associated with no response to prepulse


temp = load(elecRespPath);
elecResp = temp.elecResp;


nMovies = length(elecResp.stimInfo.movieNos);
subsetSuccessRates = ones(nMovies, 1)*inf;


elecRespFake.stimInfo = elecResp.stimInfo;

elecRespFake.analysis.latencies = cell(nMovies, 1);
elecRespFake.analysis.type = cell(nMovies, 1);
latencies = cell(nMovies, 1);

for i = 1:nMovies
    if ~isempty(subsetBin{i}) && ~isempty(elecResp.analysis.type{i})
        if sum(subsetBin{i}) > 4 %otherwise not enough to get good estimate of actual response prob
            subsetSuccessRates(i) = sum(elecResp.analysis.latencies{i}(subsetBin{i})~=0)/sum(subsetBin{i});
            elecRespFake.analysis.latencies{i} = elecResp.analysis.latencies{i}(subsetBin{i});
            elecRespFake.analysis.type{i} = 'linkage-based'; %flags movie for use in bootstrapping
        else
        end
        latencies{i} = elecResp.analysis.latencies{i}(subsetBin{i});
    end
end

% fit erf to data--subset without spike to prepulse
data = zeros(2, nMovies);
data(1,:) = elecResp.stimInfo.stimAmps;
data(2,:) = subsetSuccessRates;
lockedAmps = elecResp.analysis.finalized;

for i = nMovies: -1: 1
    if isinf(subsetSuccessRates(i))
        data(:,i) = [];
        lockedAmps(i) = [];
    end
end

data(1,:) = abs(data(1,:));
erfParams = erfFitter(data, 2, -1, 'makePlot', 0, 'lockedAmps', lockedAmps);
threshold = -erfParams(2)/erfParams(1);


thresholdStd = bootstrapThresh(elecRespFake, 100, 'scaleType', 'linear');

