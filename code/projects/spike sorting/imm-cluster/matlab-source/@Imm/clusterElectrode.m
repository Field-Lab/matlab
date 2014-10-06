function [imm, mapIndex] = clusterElectrode(imm, e)
%CLUSTERELECTRODE Run Imm algorithm on electrode e
%   imm = IMM(params, data) Generate Markov Chain and sample from it
%   for a predefined number of iterations. Save resultant set of clustering
%   solutions to the imm object.
%
%   tamachado@salk.edu 1/28/08

path       = imm.params.crpPath;
nPoints    = imm.params.nPoints;
iterations = imm.params.nIterations;
DIMS       = imm.params.nDimensions;
a_0        = imm.params.a_0;
b_0        = imm.params.b_0;
mu_0       = imm.params.mu_0;
lambda_0   = imm.params.lambda_0;
k_0        = imm.params.k_0;
v_0        = imm.params.v_0;

fprintf('Starting IMM clustering on %d points in %d dimensions...\n', nPoints, max(DIMS))

% Get projections from electrode e
d = toStruct(imm.data, e);
projections = d.projections';
nSpikes     = d.spikeCount;
clear d;

% Change to path that contains crp code
cwd = pwd;
cd(path);

% Reset BOTH random number generators
seed = 634;
rand ('state', seed);
randn('state', seed);

% Acquire random datapoints from the projections
randPoints = unidrnd(nSpikes, 1, nPoints);

% Assign the selected random points
sample = zeros(nPoints, length(DIMS));

for point = 1:nPoints
    sample(point,:) = projections(randPoints(point),DIMS);
end

try
    % Create string to write diagnostic figure to
    oPath = [imm.params.oPath '/debug_plots/' imm.params.itName sprintf('-%d', e) '.pdf'];
    
    % Check if burn in period is valid
    if(imm.params.nIterations <= imm.params.nBurnIn)
        disp('Warning: Burn in range exceeds nIterations. This may indicate a malformed config file!')
        imm.params.nBurnIn = 1;
    end

    % Run the sampler
    [classId, meanRecord, covRecord, kRecord, lpRecord, alphaRecord, bestIndex, bestIndices] = sampler(sample', iterations, a_0, b_0, mu_0, k_0, v_0, lambda_0, imm.params.figures, imm.params.nBest, imm.params.verbose, imm.params.nBurnIn, oPath);

catch
    disp('MCMC Clustering Failed!')
    cd(sprintf(cwd));
    rethrow(lasterror)
    return;
end

% Save out information to imm object
for k = 1:imm.params.nBest
    imm.clusters{e}.nClusters(k) = size(covRecord{bestIndices(k)}, 3);
end

imm.clusters{e}.classId = classId;
imm.clusters{e}.meanRecord = meanRecord;
imm.clusters{e}.covRecord = covRecord;
imm.clusters{e}.kRecord = kRecord;
imm.clusters{e}.lpRecord = lpRecord;
imm.clusters{e}.alphaRecord = alphaRecord;
imm.clusters{e}.mapIndex = bestIndex;
imm.clusters{e}.bestIndices = bestIndices;

% BestIndex contains the index of the maximum a posteriori solution
% (estimate), while bestIndices returns k number of indices
% corresponding to the k most probable solutions with unique numbers of
% clusters

% Return mapIndex
mapIndex = imm.clusters{e}.mapIndex;

% Electrode e has completed clustering successfully
imm.error(e) = imm.SUCCESS_CLUSTER;

% Restore the old working directory
cd(sprintf(cwd));