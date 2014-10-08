function elecRespArray = fitErfsConstrainedSlopeWrapper(patternNos, neuronID, pathToData, nBootstrapReps, excludePatterns)

currentDir = pwd;

cd(pathToData)

nPatterns = length(patternNos);

%% getting necessary info

analyzedMoviesAll = cell(nPatterns, 1);
allData           = cell(nPatterns-length(excludePatterns), 1);
nPulsesAll        = cell(nPatterns, 1);
nMovies           = zeros(nPatterns, 1);

pIndex = 0;
patternNosUsed = [];
for i = 1:nPatterns
    
    temp = load(['elecResp_n' num2str(neuronID) '_p' num2str(patternNos(i))]);
    elecResp = temp.elecResp;

    nMovies(i) = length(elecResp.stimInfo.movieNos);

    analyzedMoviesAll{i} = ones(nMovies(i), 1);
    for j = 1:nMovies(i)
        if isempty(elecResp.analysis.type{j})
            analyzedMoviesAll{i}(j) = 0;
        end
    end
    
    for j = 1:nMovies(i)
        nPulsesAll{i}(j) = length(elecResp.analysis.latencies{j});
    end

    if ~any(excludePatterns == patternNos(i))
        pIndex = pIndex+1;
        % stores success rate data for simultaneous erf fitting
        allData{pIndex} = zeros(2, nMovies(i));
        allData{pIndex}(1,:) = log(abs(elecResp.stimInfo.stimAmps));
        allData{pIndex}(2,:) = elecResp.analysis.successRates;
        for j = length(elecResp.stimInfo.stimAmps): -1: 1
            if isempty(elecResp.analysis.type{j})
                allData{pIndex}(:,j) = [];
            end
        end
        patternNosUsed = [patternNosUsed patternNos(i)]; %#ok<AGROW>

    end
end




%% constraining all erf slopes to be identical

[erfParamsFixedSlope threshFixedSlope errorsFixedSlope] = fitErfsWithFixedSlope(allData);

elecRespArray = cell(length(patternNosUsed), 1);
pIndex = 0;
for i = 1:nPatterns
    if ~any(excludePatterns == patternNos(i))
        pIndex = pIndex + 1;
        temp = load(['elecResp_n' num2str(neuronID) '_p' num2str(patternNos(i))]);
        elecResp = temp.elecResp;
        
        elecResp.analysis.constrainedSlopeLogThresh = threshFixedSlope(pIndex);
        elecResp.analysis.details.constrainedSlopeLogThreshPatternsUsed = patternNosUsed;
        elecResp.analysis.constrainedSlopeLogErfErr = errorsFixedSlope(pIndex);
        elecResp.analysis.constrainedSlopeLogErfParams = erfParamsFixedSlope{pIndex};

        elecRespArray{pIndex} = elecResp;
    end
end


%% bootstrap to determine standard deviations of thresholds in constrained erf slope case

newAllData = cell(length(patternNosUsed), 1);
newThreshes = cell(length(patternNosUsed), 1);
for i = 1:length(patternNosUsed)
    newThreshes{i} = zeros(nBootstrapReps, 1);
end

for j = 1:nBootstrapReps
    disp(['bootstrap repetition ' num2str(j)])
    
    %generates new set of data by sampling with replacement at each pattern/movie number to
    %determine new success rate values
    pIndex = 0;
    for k = 1:nPatterns
        if ~any(excludePatterns == patternNos(k))
            pIndex = pIndex+1;
            
            
            dataIndex = 0;
            newAllData{pIndex} = zeros(size(allData{pIndex}));
            newAllData{pIndex}(1,:) = allData{pIndex}(1,:);
            for i = 1:nMovies(k)
                if analyzedMoviesAll{k}(i)
                    dataIndex = dataIndex+1; %because allData only includes analyzed movies
                    successRate = allData{pIndex}(2, dataIndex); %measured success rate

                    nPulses = nPulsesAll{k}(i);
                    randPick = rand([nPulses, 1]);

                    while any(randPick==successRate)
                        sameInd = find(randPick==successRate);
                        for m = 1:length(sameInd);
                            randPick(sameInd(m)) = rand(1);
                        end
                    end
                    newSuccesses = sum(randPick < successRate);
                    newAllData{pIndex}(2, dataIndex) = newSuccesses/nPulses;
                end
            end
        end
    end

    % determines threshold values based on resampled data
    [erfParams newThreshesTemp] = fitErfsWithFixedSlope(newAllData);
    for k = 1:length(patternNosUsed)
        newThreshes{k}(j) = newThreshesTemp(k);
    end
end

save('bootStrap_threshholds_locked_slope.mat', 'newThreshes')

pIndex = 0;
for i = 1:nPatterns
    if ~any(excludePatterns == patternNos(i))
        pIndex = pIndex + 1;
        elecRespArray{pIndex}.analysis.constrainedSlopeLogThreshStd = std(newThreshes{pIndex});
        elecRespArray{pIndex}.analysis.details.constrainedSlopeBootstrapReps = nBootstrapReps;
    end
end

%%

cd(currentDir)