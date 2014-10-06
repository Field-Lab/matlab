function artifact = estArtWithLinkage(data, linkWindow, linkThresh, nBranches, templateFull, shiftStart, shiftEnd, shiftStep, centerChannelIndex)

% data should be 3-dimensional: traces x electrodes x samples
%
% not yet usable with multiple templates
%
% if there are multiple good electrodes, estArtWithLinkage chooses which traces to generate
% estimated artifact with separately (not necessary the same set of pulses used for each electrode's
% estimated artifact)



%linkData = data(:, linkWindow(1):linkWindow(2));

nElecs = size(data, 2);
nSamples = size(data, 3);
artifact = zeros(nElecs, nSamples);

for m = 1:nElecs
    linkData = squeeze(data(:, m, linkWindow(1):linkWindow(2)));
    template{1} = templateFull{1}(m, :); %#ok<AGROW>
    
    dist = pdist(linkData);
    link = linkage(dist, 'single');

    branches = findSubThreshBranches(link, linkThresh, nBranches);


    %finds average of each branch and number of actual branches found below threshold
    meanBranches = zeros(nBranches, nSamples);

    nBranchesFound = nBranches;
    for i = 1:nBranches
        if ~isempty(branches{i})
            %meanBranches(i, :) = mean(data(branches{i}, :), 1);
            meanBranches(i, :) = mean(data(branches{i}, m, :), 1);
        else
            nBranchesFound = i-1;
            meanBranches = meanBranches(1:nBranchesFound,:);
            break
        end
    end


    % determines which of the branches is most likely to comparing the difference between the branches
    % to the template offset at various latencies

    if nBranchesFound > 1
        lowestError = inf;
        for i = 1:nBranchesFound
            for j = 1:nBranchesFound
                testTrace = squeeze(meanBranches(i,:)-meanBranches(j,:));
                subtractedTraces = subtractWithShifts(testTrace, template, shiftStart, shiftEnd, shiftStep);

                errors = zeros(1, length(subtractedTraces));
                for k = 2:length(subtractedTraces) %because first subtracted trace corresponds to no subtraction
                    errors(k-1) = norm(subtractedTraces{k}(linkWindow(1):linkWindow(2)));
                end
                if min(errors) < lowestError
                    lowestError = min(errors);
                    bestArtifact = j;
                end
            end
        end

        artifact(m,:) = meanBranches(bestArtifact, :);

    elseif nBranchesFound == 1
        artifact(m,:) = meanBranches(1, :);
    elseif nBranchesFound == 0
        warning('No branches found below threshold: using mean of all traces as artifact estimate')
        artifact(m,:)  =  mean(data(:, m, :), 1);
    end
end
    