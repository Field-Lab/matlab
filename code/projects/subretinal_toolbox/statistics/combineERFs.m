function combineERFs(outputFolderStatistics,allSourceFolders)
% This function combines the eRFs movies computed for different pixel 
% positions into a single eRF movie structure, and plots the result.
%
% Parameters:
%   - outputFolderStatistics: statistics will be saved here
%   - allSourceFolders: a cell with all source statistics folders
%
% Returns:
%   []
%

% Version: 1.0 - 09/09/2013
% Author: Georges Goetz - Stanford University.

if nargin<2
    error('Matlab:PlotERFs:combineAndPlotERFs','Please provide at least one source statistics folder');
end
% allSourceFolders = allSourceFolders{1};
nSourceFolders = length(allSourceFolders);


%% Formatting the folders 

if outputFolderStatistics(end:end)~=filesep
    outputFolderStatistics = [outputFolderStatistics filesep];
end
for kk=1:nSourceFolders
    if allSourceFolders{kk}(end:end)~=filesep
        allSourceFolders{kk} = [allSourceFolders{kk} filesep];
    end
end

% If output folders don't exist, make them here
if ~exist(outputFolderStatistics,'dir')
    mkdir(outputFolderStatistics);
end

%% Building the neuronlist

neuronList = [];
for kk=1:nSourceFolders
    neuronList = unique(union(neuronList, readNeuronNames(allSourceFolders{kk})));
end
nNeurons = length(neuronList);

%% Combine the neurons eRFs

for kk=1:nNeurons
    neuronERFdata = combineNeuronERFs(neuronList(kk), allSourceFolders);

    % Save the eRF
    save(fullfile(outputFolderStatistics, sprintf('neuron%d_erf',neuronList(kk))),'neuronERFdata');
    
end

end % combineERFMovies

function neuronList = readNeuronNames(statsFolder)
% This function returns a list of all the neuron IDs found in a statistics
% folder.

contentsActivationStatsFolder = dir(statsFolder);
neuronList = [];

for kk=1:length(contentsActivationStatsFolder)
    if strfind(contentsActivationStatsFolder(kk).name,'.mat')
        neuronList(end+1) = str2double(contentsActivationStatsFolder(kk).name(7:end-8));
    end
end

end

function neuronERFdatacombined = combineNeuronERFs(neuronID, allSourceFolders)
% This function combines all the erf data into a single structure.

% Get an existing structure from one of the source folders and load it.
for kk=1:length(allSourceFolders)
    try load(fullfile(allSourceFolders{kk},sprintf('neuron%d_erf',neuronID)));
        break;
    catch %#ok<CTCH>
    end
end
neuronERFdatacombined = neuronERFdata;

% Now that we've got this stucture, fill it in
for kk=1:length(allSourceFolders)
    try load(fullfile(allSourceFolders{kk},sprintf('neuron%d_erf',neuronID)));
        neuronERFdatacombined.allPSTHs(neuronERFdata.pixelIndex,:) = ...
            neuronERFdata.allPSTHs(neuronERFdata.pixelIndex,:);
        neuronERFdatacombined.allResponses(neuronERFdata.pixelIndex) = ...
            neuronERFdata.allResponses(neuronERFdata.pixelIndex);
        neuronERFdatacombined.pixelIndex = union(neuronERFdatacombined.pixelIndex,...
            neuronERFdata.pixelIndex);
    catch %#ok<CTCH>
    end
end

% Fits need to be computed again: they are now incorrect.
if sum(neuronERFdatacombined.allResponses)
    [neuronERFdatacombined.mu neuronERFdatacombined.sigma] = fitGaussianToERF(neuronERFdatacombined);
else
    neuronERFdatacombined.mu = NaN;
    neuronERFdatacombined.sigma = NaN;
end

end % combineNeuronERFs

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
fitoptions = optimset('TolFun',1e-8,'TolX',1e-8);

% Solving the non-linear least-squares problem
initPos = [a b c x0 y0];
try
    [xres, resnorm, residual, exitflag] = ...
        lsqcurvefit(@fitGauss2d, initPos, [xPos yPos], zPos, ...
        [0 0 0 min(xPos) min(yPos)], [100 100 100 max(xPos) max(yPos)], fitoptions);

    mu = [xres(4) xres(5)];
    sigma = [xres(1) xres(2); xres(2) xres(3)];
catch FitError
    mu = NaN*[1 1];
    sigma = NaN*[1 1; 1 1];
    
    warning(FitError.identifier,FitError.message);
end

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