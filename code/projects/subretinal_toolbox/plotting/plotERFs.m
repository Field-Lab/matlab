function plotERFs(eRFStatsFolder, outputFolderFigures, varargin)
% plotERFMovie(eRFStatsFolder, outputFolderFigures);
%
% This function plots the static electric receptive 
% field data for all neurons found in a preparation. 
% The static ERFs are saved as an epsc file, unless specified otherwise.
%
% Parameters:
%   - eRFStatsFolder: folder in which the activation data for all the
%   neurons is saved. Activation data for a neuron with ID kkk should be
%   saved in file neuronkkk_act.mat.
%   - outputFolderFigures: folder in which the RF data figures are saved.
%
% Optional value/parameter pairs:
%   - imageFormat: format in which the plots will be saved. By default
%   'epsc'.
%
% Returns:
%   []
%

% Version: 1.0 - 09/09/2013
% Author: Georges Goetz - Stanford University.

%% Read optional inputs, format the imput folders

imageFormat = 'epsc';

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
        case 'imageformat'
            imageFormat = varargin{kk*2};
        otherwise
            err = MException('MATLAB:InvArgIn',...
                'Unknown parameter specified');
            throw(err);
    end
end

% Formatting the folders
if eRFStatsFolder(end:end)~=filesep
    eRFStatsFolder = [eRFStatsFolder filesep];
end
if outputFolderFigures(end:end)~=filesep
    outputFolderFigures = [outputFolderFigures filesep];
end

% If output folder doesn't exist, make them here
if ~exist(outputFolderFigures,'dir')
    mkdir(outputFolderFigures);
end

%% Finding all the neuron data

contentsERFStatsFolder = dir(eRFStatsFolder);
neuronNames = struct('name','');
nNeurons = 0;

for kk=1:length(contentsERFStatsFolder)
    if strfind(contentsERFStatsFolder(kk).name,'.mat')
        nNeurons = nNeurons + 1;
        neuronNames(nNeurons).name = contentsERFStatsFolder(kk).name;
    end
end

%% For each of those neurons, compute the RF, fit a gaussian, plot and save.

fh = figure();

for kk=1:nNeurons
    % Load the ERF data
    load(fullfile(eRFStatsFolder, neuronNames(kk).name));
    cNeuronID = str2double(neuronNames(kk).name(7:end-8));

    plotSingleNeuronReceptiveField(neuronERFdata, fh, cNeuronID);
    saveas(fh,fullfile(outputFolderFigures, neuronNames(kk).name(1:end-8)),imageFormat);
end

end % computeElectricReceptiveField

function plotSingleNeuronReceptiveField(neuronERFdata, fh, neuronID)
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
%     scatter(neuronERFdata.mu(1), neuronERFdata.mu(2), 'r', 'MarkerFaceColor',[1 0 0]);
end
% If the neuron position was specified, add a marker
if ~isempty(neuronERFdata.neuronPos)
    scatter(neuronERFdata.neuronPos(1), neuronERFdata.neuronPos(2), 'b', 'MarkerFaceColor',[0 0 1]);
end
% Title and stuff
axis equal; axis off;
title(sprintf('Neuron %d electric receptive field; max #spikes elicited: %2.2f', ...
    neuronID, max(neuronERFdata.allResponses)),...
    'fontsize',14)

end % plotSingleNeuronReceptiveField

function drawEllipse(mu,V,color)
% Draws a simple ellipse

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

end % drawEllipse