function plotSlopeVArraySpace(cellInfo, varargin)


% cellInfo should be a struct array with the following fields:
%
%	.pathToAxonDirData = '/snle/lab/Experiments/Array/Analysis/2008-08-26-0/data007-NW/axonDirData.mat';
%   .details.lambdas
%   .details.sElecs %with order corresponding to lambdas


p = inputParser;

p.addRequired('cellInfo', @isstruct)

p.addParamValue('excludeBadElecs', true, @islogical)
p.addParamValue('binSpacing', 0.1, @isnumeric) %whether to normalize lambdas between 0 and 1 for each cell individually or for all cells together

p.parse(cellInfo, varargin{:})

binSpacing = p.Results.binSpacing;
excludeBadElecs = p.Results.excludeBadElecs;

allLam60 = [];
allLam30 = [];
for ii = 1:length(cellInfo)
    if cellInfo(ii).pitch == 60;
        if excludeBadElecs
            allLam60 = [allLam60; cellInfo(ii).details.lambdas(~ismember(cellInfo(ii).details.sElecs, cellInfo(ii).badSElecs))'];
        else
            allLam60 = [allLam60; cellInfo(ii).details.lambdas];

        end
    else
        if excludeBadElecs
            allLam30 = [allLam30; cellInfo(ii).details.lambdas(~ismember(cellInfo(ii).details.sElecs, cellInfo(ii).badSElecs))'];
        else
            allLam60 = [allLam60; cellInfo(ii).details.lambdas];
        end
    end
end

maxLam = max([max(allLam60) max(allLam30)]);
minLam = min([min(allLam60) min(allLam30)]);

minLam = binSpacing*(floor(minLam*(1/binSpacing))-1); %round down to nearest binSpacing multiple - 1
maxLam = binSpacing*(ceil(maxLam*(1/binSpacing))+1); %round up to nearest binSpacing multiple + 1

binEdges = minLam:binSpacing:maxLam;

histVals60 = histc(allLam60, binEdges); %hist val of last bin should be 0
histVals30 = histc(allLam30, binEdges); %hist val of last bin should be 0

%plot histogram
patchVals30 = [binEdges(1); 0];
patchVals60 = [binEdges(1); 0];

for ii = 1:length(binEdges)-1
    patchVals30 = [patchVals30, [binEdges(ii) binEdges(ii+1); histVals30(ii) histVals30(ii)]];
    patchVals60 = [patchVals60, [binEdges(ii) binEdges(ii+1); histVals60(ii) histVals60(ii)]];
end
patchVals30 = [patchVals30, [binEdges(end); 0]];
patchVals60 = [patchVals60, [binEdges(end); 0]];


figure
hold on
patch(patchVals30(1,:), patchVals30(2,:)/length(allLam30), [1 0 0])
patch(patchVals60(1,:), patchVals60(2,:)/length(allLam60), [0 0 1])

