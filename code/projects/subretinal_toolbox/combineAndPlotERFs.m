function combineAndPlotERFs(outputFolderFigures,outputFolderStatistics,allSourceFolders,...
                               neuronPos,varargin)
% This function combines the eRFs computed for different pixel positions
% into a single eRF structure, and plots the result.
%
% Parameters:
%   - outputFolderFigures: figures will be saved here
%   - outputFolderStatistics: statistics will be saved here
%   - allSourceFolders: a cell with all source statistics folders
%   - neuronPos: optional, position of the neurons on the array.
%
% Returns:
%   []
%

if nargin<3
    error('Matlab:PlotERFs:combineAndPlotERFs','Please provide at least one source statistics folder');
end
% allSourceFolders = allSourceFolders{1};
nSourceFolders = length(allSourceFolders);
imageFormat = 'epsc';
fh = figure();
if ~exist('neuronPos','var')
    neuronPos = [];
end

%% Formatting the folders 

if outputFolderStatistics(end:end)~=filesep
    outputFolderStatistics = [outputFolderStatistics filesep];
end
if outputFolderFigures(end:end)~=filesep
    outputFolderFigures = [outputFolderFigures filesep];
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
if ~exist(outputFolderFigures,'dir')
    mkdir(outputFolderFigures);
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
    
    % If the neuron position was specified, find it
    if ~isempty(neuronPos)
        cNeuronPos = neuronPos(kk, :);
    else
        cNeuronPos = [];
    end

    % Save the eRF
    save(fullfile(outputFolderStatistics, sprintf('neuron%d_erf',neuronList(kk))),'neuronERFdata');
    
    % Plot and save
    plotSingleNeuronReceptiveField(neuronERFdata, fh, neuronList(kk), cNeuronPos);
    saveas(fh,fullfile(outputFolderFigures, sprintf('neuron%d',neuronList(kk))),imageFormat);
    
end

end % combineAndPlotERFs

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
        neuronERFdatacombined.allResponses(neuronERFdata.pixelIndex) = ...
            neuronERFdata.allResponses(neuronERFdata.pixelIndex);
        neuronERFdatacombined.pixelIndex = union(neuronERFdatacombined.pixelIndex,...
            neuronERFdata.pixelIndex);
    catch %#ok<CTCH>
    end
end

if sum(neuronERFdatacombined.allResponses)
    [neuronERFdatacombined.mu neuronERFdatacombined.sigma] = fitGaussianToERF(neuronERFdatacombined);
else
    neuronERFdatacombined.mu = NaN;
    neuronERFdatacombined.sigma = NaN;
end

end

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