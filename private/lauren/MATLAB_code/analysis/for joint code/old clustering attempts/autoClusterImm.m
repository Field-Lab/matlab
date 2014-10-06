function clusterIndex = autoClusterImm(projections, verbose, iterations, DIMS, nBurnIn, a_0, b_0, mu_0, k_0, lambda_0, v_0)
% adapted from clusterElectrode by Tim Machado (1/28/08)
%
%AUTOCLUSTERIMM Runs Imm algorithm on data in projections
%   Generates Markov Chain and samples from it
%   for a predefined number of iterations.
%
%   Contains all parameters that control how the clustering algorithm is run.
%
% arguments
%   projections: an array with points x dimensions
%
%
%


% path to sampler.m
path = '/snle/lab/Development/RRS/imm-cluster/crp-source/crp';
addpath '/snle/lab/Development/RRS/imm-cluster/crp-source/crp' '/snle/lab/Development/RRS/imm-cluster/crp-source/utilities'

% parameters (values from config.m)
figures = false; %specifies whether figures of cluster number, probability, and value of alpha vs. iteration are shown
writeFigure = false; %specifies whether these figures are written to disk
itName = '2008-08-26-0-data000'; %name of figure saved
oPathPrefix = ['/Volumes/Creampuff/Analysis.noindex/Lauren/2008-08-26-0/data000-imm/']; %location where figures will be written to

nBest = 1; %specifies the number of maximum probability solutions (with distinct numbers of clusters) that should be found

if nargin < 2
    verbose = true; %specifies whether test output is displayed during clustering
elseif nargin < 3
    % Iteration settings (values from config.m)
    iterations =  250;
elseif nargin < 4
    DIMS = 5; %number of PCs included in clustering
elseif nargin < 5
    nBurnin = 125;
% CRP/IMM Hyperparameters (values from config.m)
% Use sample_igmm_prior.m to determine the best priors to use
% Parameters for the gamma function (sampled from when computing 
% alpha in the chinese restaurant process sampler)
elseif nargin < 6
    a_0 = 1; %increasing this value increases the number of clusters predicted by priors
elseif nargin < 7
    b_0 = 1; %increasing this value increases the number of clusters predicted by priors
% Priors on cluster means
elseif nargin < 8
    mu_0 = zeros(DIMS,1);
elseif nargin < 9
    k_0 = 0.1;
% Priors on cluster covariances
elseif nargin < 10
    lambda_0 = eye(DIMS)*1000;
elseif nargin < 11
    v_0 = 10;
end

%%
nSpikes = size(projections, 2);

fprintf('Starting IMM clustering on %d points in %d dimensions...\n', nSpikes, DIMS)

% Change to path that contains crp code
cwd = pwd;
cd(path);

% Reset BOTH random number generators
seed = 634;
rand ('state', seed);
randn('state', seed);

% take specified number of dimensions from data
sample = projections(:, 1:DIMS);

% try
    
    % Check if burn in period is valid
    if(iterations <= nBurnIn)
        disp('Warning: Burn in range exceeds nIterations. This may indicate a malformed config file!')
        nBurnin = 1;
    end

    % Run the sampler
    
    if writeFigure
        % Create string to write diagnostic figure to
        oPath = [oPathPrefix '/debug_plots/' itName '.pdf'];
        % run clustering algorithm
        [classId, meanRecord, covRecord, kRecord, lpRecord, alphaRecord, bestIndex, bestIndices] = sampler(sample', iterations, a_0, b_0, mu_0, k_0, v_0, lambda_0, figures, nBest, verbose, nBurnIn, oPath);
    else 
        [classId, meanRecord, covRecord, kRecord, lpRecord, alphaRecord, bestIndex, bestIndices] = sampler(sample', iterations, a_0, b_0, mu_0, k_0, v_0, lambda_0, figures, nBest, verbose, nBurnIn);
    end
        
% catch
%     disp('MCMC Clustering Failed!')
%     cd(sprintf(cwd));
%     rethrow(lasterror)
%     return;
% end

% classId: record of cluster assignments
% meanRecord: record of cluster means at each iteration
% covRecord: record of cluster covariances at each iteration
% kRecord: number of clusters at each iteration
% lpRecord: log probability
% bestIndex: index of the iteration that has maximum log probability
% bestIndices: if nBest>1, gives nBest indeces of iterations that have maximum log probability and
% unique numbers of clusters

clusterIndex = classId(:,bestIndex);


% Restore the old working directory
cd(sprintf(cwd));