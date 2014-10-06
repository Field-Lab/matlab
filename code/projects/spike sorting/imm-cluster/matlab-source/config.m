function params = config(range, itName, prjPath, oPath, modelPath, nProcess, params, nElectrodes, noOutput)
%CONFIG Script to generate config.mat with specified parameters
%   User should specify the four arguments to this function:
%      1. range     -- Range of electrodes to cluster e.g. 1:512
%      2. name      -- Name of this run e.g. 2005-02-06-1-data010-part1
%      3. prjPath   -- Path to projections file to cluster
%      4. oPath     -- Path to directory in which to store output
%      5. modelFile -- Path to model file containing eigenvectors
%
%   Other parameters that are varied infrequently can be modified
%   by actually editing config.m
%
%   After running this, simply type ImmCluster to begin clustering
%
%   tamachado@salk.edu 2/5/08

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check if a params struct was passed in
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist('params', 'var') && isa(params, 'struct')
    existParams = true;
else
    existParams = false;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters to change each run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Iteration name (used as a prefix for all output files)
if nargin < 2
    params.itName = '2007-02-06-1-data015';
else
    params.itName = itName;
end

% Path to .prj file for dataset to analyze
if nargin < 3
    params.prjPath = '/Volumes/Sleepy/Data/Machado/2007-02-06-1/data015/data015.prj';
else
    params.prjPath = prjPath;
end

% Path to output directory
if nargin < 4
    params.oPath = ['/Volumes/Sleepy/Data/Machado/matlab/output/ce/']; %params.itName];
else
    params.oPath = oPath;
end

% Range of electrodes to cluster
if nargin < 1
    params.range = 1:1;
else
    params.range = range;
end

% Number of MATLAB processes used to do IMM clustering
if nargin < 6
    params.nProcess = 4;
else
    params.nProcess = nProcess;
end

% Number of electrodes in this dataset
if nargin < 8
    params.nElectrodes  = 512;
else
    params.nElectrodes = nElectrodes;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters that rarely change
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Path to immUtilities and crp code
javaaddpath /snle/lab/Development/RRS/vision-package/vision.app/Contents/Resources/Java/Vision.jar;
addpath('/snle/lab/Development/RRS/imm-cluster/crp-source/utilities');
addpath('/snle/lab/Development/RRS/imm-cluster/crp-source/distributions');

% Path to .model file containing eigenvectors
if nargin < 5 || ~isa(modelPath, 'char')
    params.modelPath = 'false';
    params.modelOnly = 'false';
else
    params.modelOnly = 'true';
    params.modelPath = modelPath;
end

% Path to CRP/IMM sampler code
params.crpPath = '/snle/lab/Development/RRS/imm-cluster/crp-source/crp';

% Do clusters already exist in clusters.mat?
params.clustersExist = false;

% Only set default parameters if none of them are set already
if existParams == false
    
    % Iteration settings
    params.nIterations  = 250;
    params.nDimensions  = 1:5;
    params.nBurnIn      = 125;
    params.nPoints      = 1200;
    params.nBest        = 1; %note this!

    % CRP/IMM Hyperparameters
    % Use sample_igmm_prior.m to determine the best priors to use
    params.a_0      = 1;
    params.b_0      = 1;
    params.mu_0     = zeros(length(params.nDimensions),1);
    params.lambda_0 = eye(length(params.nDimensions))* 1000;
    params.k_0      = .1;
    params.v_0      = 10;

    % Smallest cluster to add to neurons file
    params.minSpikes = 5;

    % Should the IMM algorithm display figures?
    params.figures   = true;

    % Should the IMM algorithm display any text output?
    params.verbose   = true;

    % Should this script write plots to the output directory?
    params.savePlots = false;
end

% Make necessary directories
mkdir([params.oPath]);
mkdir([params.oPath '/neurons/']);
mkdir([params.oPath '/model/']);
mkdir([params.oPath '/debug_plots/']);
mkdir([params.oPath '/plots/']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save config.mat file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If noOutput is specified, do not actually write config to a file
if nargin < 9
    path = 'config.mat';
    save( path, 'params');
end
