function [erfParams thresholds errors] = fitErfsWithFixedSlope(allData, varargin)

%allData{i}(1,:) should be log(stimulus amplitudes)


p = inputParser;

p.addRequired('allData', @(x)isa(x, 'cell'))

p.addParamValue('iterateFits', 0, @(x)any(x==[0 1]))



p.parse(allData, varargin{:})

iterateFits = p.Results.iterateFits;


%%

errordlg('This function needs to be updated to reflect change in erfFitter to using maximum likelihood')



nPatterns = length(allData);

adjustedData = cell(nPatterns, 1);

%figure
%hold on
catData = [];
for i = 1:nPatterns
    erfParamsInit = erfFitter(allData{i}, 2, 1, 'makePlot', 0);
    
    logThresh = -erfParamsInit(2)/erfParamsInit(1);
    
    %xProj = log(0.1:0.001:1);
    %projection = 0.5 + 0.5*erf(erfParamsInit(1)*xProj+erfParamsInit(2));
    
    %plot(xProj-logThresh, projection)
    
    adjustedData{i} = allData{i};
    adjustedData{i}(1,:) = allData{i}(1,:) - logThresh;
    %plot(adjustedData{i}(1,:), adjustedData{i}(2,:))
    %plot(allData{i}(1,:), allData{i}(2,:), 'r')
    catData = cat(2, catData, adjustedData{i});
end
%hold off

erfParamsAll = erfFitter(catData, 2, 1, 'makePlot', 1);


erfParams = cell(nPatterns, 1);
thresholds = zeros(nPatterns, 1);
errors = zeros(nPatterns, 1);

for i = 1:nPatterns
    [erfParams{i} x errors(i)] = erfFitter(allData{i}, erfParamsAll(1), 1, 'makePlot', 0, 'setParam', 1);
    thresholds(i) = -erfParams{i}(2)/erfParams{i}(1);
end

if iterateFits
    difference = 100;
    iterationCounter = 0;
    iteratedThresholds = thresholds;
    while difference > 1e-5
        iterationCounter = iterationCounter+1
        adjustedData = cell(nPatterns, 1);
        for i = 1:nPatterns
            adjustedData{i} = allData{i};
            adjustedData{i}(1,:) = allData{i}(1,:) - iteratedThresholds(i, iterationCounter);
            catData = cat(2, catData, adjustedData{i});
        end
        erfParamsAll = erfFitter(catData, 2, 1, 'makePlot', 1);
        
        for i = 1:nPatterns
            erfParamsTemp = erfFitter(allData{i}, erfParamsAll(1), 1, 'makePlot', 0, 'setParam', 1);
            iteratedThresholds(i, iterationCounter+1) = -erfParamsTemp(2)/erfParamsTemp(1);
        end
        
        difference = norm(squeeze(iteratedThresholds(:, iterationCounter+1))-squeeze(iteratedThresholds(:, iterationCounter)))
        
    end
end
