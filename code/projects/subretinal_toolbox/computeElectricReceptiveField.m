function computeElectricReceptiveField(activationStatsFolder, deviceType, pixelIndex, ...
    outputFolderStatistics, outputFolderFigures, varargin)
% computeElectricReceptiveField(activationFolder, deviceType, pixelIndex, ...
%    outputFolderStatistics, outputFolderFigures, ...);
%
% This function computes the electric receptive field data for all neurons
% found in a preparation. The receptive fields are saved in the output
% folder, in a structure containing for each neuron the following fields: 
%    - rf.pixelPos: a nPixels*2 matrix, in which each row represents the
%   position of the center of a pixel on the array
%    - rf.allResponses: number of spikes elicited per trial for each of the
%   pixels
%    - rf.mu: mean position of the RF
%    - rf.sigma: covariance matrix for the RF
%    - rf.surroundPixels: a structure containing position of surrounding
%   pixels, which is used only for plotting.
%
% Parameters:
%   - activationFolder: folder in which the activation data for all the
%   neurons is saved. Activation data for a neuron with ID kkk should be
%   saved in file neuronkkk_act.mat.
%   - deviceType: a string, which can be either 'small' or 'medium'. It is
%   used to link to a file which has the coordinates of the center of the
%   small or medium pixels.
%   - pixelIndex: a list of which pixels have been stimulated, so that 
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

% Version: 1.1 - 05/30/2013
% Author: Georges Goetz - Stanford University.

%% Read optional inputs, format the imput folders

% Default values
savePlots = true;
imageFormat = 'epsc';
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
        case 'saveplots'
            savePlots = varargin{kk*2};
        case 'imageformat'
            imageFormat = varargin{kk*2};
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
if activationStatsFolder(end:end)~=filesep
    activationStatsFolder = [activationStatsFolder filesep];
end
if outputFolderStatistics(end:end)~=filesep
    outputFolderStatistics = [outputFolderStatistics filesep];
end
if outputFolderFigures(end:end)~=filesep
    outputFolderFigures = [outputFolderFigures filesep];
end

% If output folders don't exist, make them here
if ~exist(outputFolderStatistics,'dir')
    mkdir(outputFolderStatistics);
end
if ~exist(outputFolderFigures,'dir')
    mkdir(outputFolderFigures);
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

contentsActivationStatsFolder = dir(activationStatsFolder);
neuronNames = struct('name','');
nNeurons = 0;

for kk=1:length(contentsActivationStatsFolder)
    if strfind(contentsActivationStatsFolder(kk).name,'.mat')
        if ~isempty(neuronList)
            cNeuronID = str2double(contentsActivationStatsFolder(kk).name(7:end-8));
            if nnz(neuronList==cNeuronID)>0
                nNeurons = nNeurons + 1;
                neuronNames(nNeurons).name = contentsActivationStatsFolder(kk).name;
            end
        else
            nNeurons = nNeurons + 1;
            neuronNames(nNeurons).name = contentsActivationStatsFolder(kk).name;
        end
    end
end

%% For each of those neurons, compute the RF, fit a gaussian, plot and save.

if savePlots
    fh = figure();
end

for kk=1:nNeurons
    % Load the activation data
    load(fullfile(activationStatsFolder, neuronNames(kk).name));
    cNeuronID = str2double(neuronNames(kk).name(7:end-8));
    
    % If the neuron position was specified, find it
    if ~isempty(neuronPos)
        cNeuronPos = neuronPos(neuronList == cNeuronID, :);
    else
        cNeuronPos = [];
    end
    
    % Compute the eRF
    neuronERFdata = computeERF(activationData, pixelIndex, pixelPos, surroundPixelPos);
    
    % Save the eRF
    save(fullfile(outputFolderStatistics, sprintf('neuron%d_erf',...
        str2double(neuronNames(kk).name(7:end-8)))),'neuronERFdata');
    
    % If saving the plot, plot and save
    if savePlots
        plotSingleNeuronReceptiveField(neuronERFdata, fh, cNeuronID, cNeuronPos);
        saveas(fh,fullfile(outputFolderFigures, neuronNames(kk).name(1:end-8)),imageFormat);
    end
end

end % computeElectricReceptiveField

function neuronERFdata = computeERF(neuronActivationData, pixelIndex, pixelPos, surroundPixelPos)
% This function computes and returns the neuron eRF data in a structure as
% described in the main function.
% 
% The fields of the structure are as follows:
%    - neuronERFdata.pixelPos: a nPixels*2 matrix, in which each row 
%   represents the position of the center of a pixel on the array
%    - neuronERFdata.allResponses: number of spikes elicited per trial for
%   each of the pixels
%    - neuronERFdata.mu: mean position of the RF
%    - neuronERFdata.sigma: covariance matrix for the RF
%    - neuronERFdata.surroundPixels: a structure containing position of 
%   surrounding pixels, which is used only for plotting.

% Getting number of elicited spikes
allResponsesMeasured = neuronActivationData.data(:,ismember(neuronActivationData.labels,'nSpikes'));
assert(numel(pixelIndex)==numel(allResponsesMeasured));

% Sorting according to the order the pixels were stimulated in
allResponses = zeros(size(pixelPos,1),1);
allResponses(pixelIndex) = allResponsesMeasured;

% Creating the result structure
neuronERFdata.pixelIndex = pixelIndex;
neuronERFdata.pixelPos = pixelPos;
neuronERFdata.allResponses = allResponses;
if sum(allResponses)
    [neuronERFdata.mu neuronERFdata.sigma] = fitGaussianToERF(neuronERFdata);
else
    neuronERFdata.mu = NaN;
    neuronERFdata.sigma = NaN;
end
neuronERFdata.surroundPixels = surroundPixelPos;

end % computeERF

function plotSingleNeuronReceptiveField(neuronERFdata, fh, neuronID, neuronPos)
% Function dedicated to plotting the eRF for a single neuron.

% Formatting...
allPixelPositions = [neuronERFdata.pixelPos; neuronERFdata.surroundPixels];
allOutputs = [neuronERFdata.allResponses; zeros(size(neuronERFdata.surroundPixels,1),1)];
if nnz(allOutputs)
    allOutputs = allOutputs/max(allOutputs);
end
nPixels = size(neuronERFdata.pixelPos,1);

[V,C] = voronoin([allPixelPositions(:,1),allPixelPositions(:,2)]);

% Plotting the result
fh = figure(fh); clf; set(fh,'color','white');
hold on
% Plotting the pixels
for ii = 1:nPixels
    if all(C{ii}~=1)   % If at least one of the indices is 1, 
                      % then it is an open region and we can't 
                      % patch that.
        patch(V(C{ii},1),V(C{ii},2),repmat(allOutputs(ii),1,3)); % use color i.
    end
end
% Plotting the RF
if ~isnan(sum(neuronERFdata.mu))
    drawEllipse(neuronERFdata.mu,neuronERFdata.sigma,'r')
    scatter(neuronERFdata.mu(1), neuronERFdata.mu(2), 'r', 'MarkerFaceColor',[1 0 0]);
end
% If the neuron position was specified, add a marker
if ~isempty(neuronPos)
    scatter(neuronPos(1), neuronPos(2), 'b', 'MarkerFaceColor',[0 0 1]);
end
% Title and stuff
axis equal; axis off;
title(sprintf('Neuron %d electric receptive field; max #spikes elicited: %2.2f', ...
    neuronID, max(neuronERFdata.allResponses)),...
    'fontsize',14)

end % plotSingleNeuronReceptiveField

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

function drawEllipse(mu,V,color)

theta = 2*pi*[0:0.01:1];
x = [cos(theta); sin(theta)];
[V,D] = eig(inv(V));
x = V*sqrt(D)*x;

for ii=1:length(mu)
    x(ii,:) = x(ii,:) + mu(ii);
end

if nargin<3
    plot(x(1,:),x(2,:));
else
    plot(x(1,:),x(2,:),color);
end

end