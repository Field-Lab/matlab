function computeERFs(psthStatsFolder, activationFolder, deviceType, pixelIndex, ...
    outputFolderStatistics, varargin)
% computeERF(psthStatsFolder, activationFolder, deviceType, pixelIndex, ...
%    outputFolderStatistics, ...);
%
% This function computes the electric receptive field data for all neurons
% found in a preparation. The receptive fields are saved in the output
% folder, in a structure containing for each neuron the following fields: 
%    - rf.pixelPos: a nPixels*2 matrix, in which each row represents the
%   position of the center of a pixel on the array
%    - rf.allResponses: PSTH for each of the pixels
%    - rf.mu: mean position of the RF
%    - rf.sigma: covariance matrix for the RF
%    - rf.surroundPixels: a structure containing position of surrounding
%   pixels, which is used only for plotting.
%
% Parameters:
%   - psthStatsFolder: folder in which the neuronInformation structures
%   have been saved.
%   - activationFolder: folder in which the activation data for all the
%   neurons is saved. Activation data for a neuron with ID kkk should be
%   saved in file neuronkkk_act.mat.
%   - deviceType: a string, which can be either 'small' or 'medium'. It is
%   used to link to a file which has the coordinates of the center of the
%   small or medium pixels.
%   - pixelIndex: a list of which pixels have been stimulated in which order. 
%   - outputFolderStatistics: folder in which the RF data statistics are
%   saved.
%   - outputFolderFigures: folder in which the RF data figures are saved.
%
% Optional value/parameter pairs:
%   - neuronList: if specified, only the IDs specified will be plotted.
%   Otherwise, every available neuronID will be plotted.
%   - neuronPos: if specified, the position of the neuron will be
%   visualised on the array. When the positions of the neurons are
%   specified, the corresponding neuronList has to be specified too,
%   otherwise there is no way of knowing which position corresponds to
%   which neuron.
%
% Returns:
%   []
%

% Version: 1.2 - 05/30/2013
% Author: Georges Goetz - Stanford University.

%% Read optional inputs, format the imput folders

% Default values
neuronList = [];
neuronPos = [];

% Checking the optional parameters
nbin = length(varargin);
if mod(nbin,2)==1
    err = MException('MATLAB:InvArgIn', ...
        'Unexpected number of arguments');
    throw(err);
end

% Reading the optional input arguments
for kk=1:(nbin/2)
    if ~ischar(varargin{kk*2-1})
        err = MException('MATLAB:InvArgIn',...
            'Unexpected additional property');
        throw(err);
    end
    
    switch lower(varargin{kk*2-1})
        case 'neuronlist'
            neuronList = varargin{kk*2};
        case 'neuronposition'
            neuronPos = varargin{kk*2};
        otherwise
            err = MException('MATLAB:InvArgIn',...
                'Unknown parameter specified');
            throw(err);
    end
end

% The neuron ID list has to be specified if the neuronPos variable has been
% specified
if isempty(neuronList)&&~isempty(neuronPos)
    error('When the neuron position is specified, the neuron ID list has to be specified too');
end

% Formatting the folders
if psthStatsFolder(end:end)~=filesep
    psthStatsFolder = [psthStatsFolder filesep];
end
if activationFolder(end:end)~=filesep
    activationFolder = [activationFolder filesep];
end
if outputFolderStatistics(end:end)~=filesep
    outputFolderStatistics = [outputFolderStatistics filesep];
end

% If output folders don't exist, make them here
if ~exist(outputFolderStatistics,'dir')
    mkdir(outputFolderStatistics);
end

%% Load the pixel position information

switch deviceType
    case 'medium'
        pixelPos = dlmread('med_array_lattice.txt');
        surroundPixelPos = dlmread('med_array_surround.txt');
    case 'small'
        pixelPos = dlmread('small_array_lattice_flipped-lr.txt');
        surroundPixelPos = dlmread('small_array_surround_flipped-lr.txt');
    case 'large'
        error('Not yet implemented');
end

%% Finding all the neuron data

contentsPsthFolder = dir(psthStatsFolder);
neuronNames = struct('name','');
nNeurons = 0;

for kk=1:length(contentsPsthFolder)
    if strfind(contentsPsthFolder(kk).name,'.mat')
        if ~isempty(neuronList)
            cNeuronID = str2double(contentsPsthFolder(kk).name(7:end-4));
            if nnz(neuronList==cNeuronID)>0
                nNeurons = nNeurons + 1;
                neuronNames(nNeurons).name = contentsPsthFolder(kk).name;
            end
        else
            nNeurons = nNeurons + 1;
            neuronNames(nNeurons).name = contentsPsthFolder(kk).name;
        end
    end
end

%% Get some basic PSTH information from at least one neuron

for kk=1:nNeurons
    try 
        load(fullfile(psthStatsFolder, neuronNames(kk).name));
        nColsPSTH = length(obj.PSTHs(1).data);
        break;
    catch %#ok<CTCH>
    end
end

if ~exist('nColsPSTH','var')
    loadException = MException('MATLAB:FileIO:InvalidFileNames',...
        sprintf('No neuron statistics found in the folder "%s"',psthStatsFolder));
    throw(loadException);
end
    

%% For each of those neurons, compute the RF, fit a gaussian, plot and save.

for kk=1:nNeurons
    cNeuronID = str2double(neuronNames(kk).name(7:end-4));
    
    % If the neuron position was specified, find it
    if ~isempty(neuronPos)
        cNeuronPos = neuronPos(neuronList == cNeuronID, :);
    else
        cNeuronPos = [];
    end
    
    % Load the activation data
    try
        load(fullfile(psthStatsFolder, neuronNames(kk).name));
        load(fullfile(activationFolder, sprintf('neuron%d_act.mat',cNeuronID)));

        % Compute the eRF
        neuronERFdata = computeSingleNeuronERF(obj, activationData, cNeuronPos, ...
            pixelIndex, pixelPos, surroundPixelPos);
        
    catch loadingError
        warning(loadingError.identifier, loadingError.message);
        
        % If loading didn't work, replace by an empty ERF.
        neuronERFdata = createEmptyTdERF(nColsPSTH, cNeuronPos, pixelIndex, pixelPos, surroundPixelPos);
    end
    
    % Save the eRF
    save(fullfile(outputFolderStatistics, sprintf('neuron%d_erf',...
            cNeuronID)),'neuronERFdata');
            
end

end % computeERFs

function neuronERFdata = computeSingleNeuronERF(neuronPsthData, neuronActivationData, ...
    neuronPos, pixelIndex, pixelPos, surroundPixelPos)
% This function computes and returns the neuron eRF data in a structure as
% described in the main function.
% 
% The fields of the structure are as follows:
%    - neuronTdERFdata.pixelPos: a nPixels*2 matrix, in which each row 
%   represents the position of the center of a pixel on the array
%    - neuronERFdata.allResponses: number of spikes elicited per trial for
%   each of the pixels
%    - neuronTdERFdata.allPSTHs: psth for each pixel activation
%    - neuronTdERFdata.mu: mean position of the RF
%    - neuronTdERFdata.sigma: covariance matrix for the RF
%    - neuronTdERFdata.surroundPixels: a structure containing position of 
%   surrounding pixels, which is used only for plotting.
%    - neuronTdERFdata.neuronPos: coordinates of the EI on the array.

% Getting the PSTHs
allPsthsMeasured = neuronPsthData.PSTHs;
assert(numel(pixelIndex)==numel(allPsthsMeasured));

% Getting the number of elicited spikes
allResponsesMeasured = neuronActivationData.data(:,ismember(neuronActivationData.labels,'nSpikes'));
assert(numel(pixelIndex)==numel(allResponsesMeasured));

% Sorting according to the order the pixels were stimulated in,  first PSTHS
allPSTHs = zeros(size(pixelPos,1),length(allPsthsMeasured(1).data));
for kk=1:length(pixelIndex)
    allPSTHs(pixelIndex(kk),:) = neuronPsthData.PSTHs(kk).data;
end
% then sorting number of spikes
allResponses = zeros(size(pixelPos,1),1);
allResponses(pixelIndex) = allResponsesMeasured;

% Creating the result structure
neuronERFdata.pixelIndex = pixelIndex;
neuronERFdata.pixelPos = pixelPos;
neuronERFdata.allPSTHs = allPSTHs;
neuronERFdata.surroundPixels = surroundPixelPos;
neuronERFdata.neuronPos = neuronPos;
neuronERFdata.allResponses = allResponses;

% If the neuron was simulated, fit a gaussian.
if sum(allResponses)
    [neuronERFdata.mu neuronERFdata.sigma] = fitGaussianToERF(neuronERFdata);
else
    neuronERFdata.mu = NaN;
    neuronERFdata.sigma = NaN;
end

end % computeSingleNeuronERF

function neuronTdERFdata = createEmptyTdERF(nColsPSTH, neuronPos, pixelIndex, pixelPos, surroundPixelPos)
% This function creates a blank Td ERF structure (= no stimulation
% anywhere), and should be called in case a neuron was not found in a
% partial ERF mapping dataset. 


allPSTHs = zeros(size(pixelPos,1),nColsPSTH);
neuronTdERFdata.pixelIndex = pixelIndex;
neuronTdERFdata.pixelPos = pixelPos;
neuronTdERFdata.allPSTHs = allPSTHs;
neuronTdERFdata.surroundPixels = surroundPixelPos;
neuronTdERFdata.neuronPos = neuronPos;

end % createEmptyTdERF

function [mu, sigma] = fitGaussianToERF(neuronERFdata)
% Fits a 2D gaussian to the eRF for a single neuron.

xPos = neuronERFdata.pixelPos(:,1);
yPos = neuronERFdata.pixelPos(:,2);
zPos = neuronERFdata.allResponses/max(neuronERFdata.allResponses);

% Starting point
x0 = sum(xPos.*zPos)/sum(zPos);
y0 = sum(yPos.*zPos)/sum(zPos);
a = 1e-10;
b = 0;
c = 1e-10;

% Setting options for the termination condition of the fit
fitoptions = optimset('TolFun',1e-6,'TolX',1e-5);

% Solving the non-linear least-squares problem
initPos = [a b c x0 y0];
[xres, resnorm, residual, exitflag] = lsqcurvefit(@fitGauss2d, initPos, [xPos yPos], zPos, ...
    [0 0 0 min(xPos) min(yPos)], [100 100 100 max(xPos) max(yPos)], fitoptions);

mu = [xres(4) xres(5)];
sigma = [xres(1) xres(2); xres(2) xres(3)];

end % fitGaussianToERF

function val = fitGauss2d(x,XY)
% Function we're trying to fit to the data
% x: vector containing the different coefficients of the equation.
% XY = [X, Y], coordinates of the points in 2D
% Z: z data set.

a = x(1);
b = x(2);
c = x(3);
x0 = x(4);
y0 = x(5);

X = XY(:,1);
Y = XY(:,2);

val = exp(-(a*(X-x0).^2 + 2*b*(X-x0).*(Y-y0) + c*(Y-y0).^2));

end % fitGauss2d