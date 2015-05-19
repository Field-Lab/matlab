function plotERFMovie(tdERFStatsFolder, outputFolderFigures, varargin)
% plotERFMovie(tdERFStatsFolder, outputFolderFigures, ...);
%
% This function plots the time-dependent electric receptive 
% field data for all neurons found in a preparation. 
% The time-dependent ERFs are saved as an animated gif.
%
% Parameters:
%   - tdERFStatsFolder: folder in which the activation data for all the
%   neurons is saved. Activation data for a neuron with ID kkk should be
%   saved in file neuronkkk_act.mat.
%   - outputFolderFigures: folder in which the RF data figures are saved.
%
% Optional value/parameter pairs:
%   - nFrames: the number of 5ms frames in the resulting eRF movie.
%
% Returns:
%   []
%

% Version: 1.0 - 08/09/2013
% Author: Georges Goetz - Stanford University.

%% Read optional inputs, format the imput folders

nFrames = 40;

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
        case 'nframes'
            nFrames = varargin{kk*2};
        otherwise
            err = MException('MATLAB:InvArgIn',...
                'Unknown parameter specified');
            throw(err);
    end
end

% Formatting the folders
if tdERFStatsFolder(end:end)~=filesep
    tdERFStatsFolder = [tdERFStatsFolder filesep];
end
if outputFolderFigures(end:end)~=filesep
    outputFolderFigures = [outputFolderFigures filesep];
end

% If output folder doesn't exist, make them here
if ~exist(outputFolderFigures,'dir')
    mkdir(outputFolderFigures);
end

%% Finding all the neuron data

contentsTdERFStatsFolder = dir(tdERFStatsFolder);
neuronNames = struct('name','');
nNeurons = 0;

for kk=1:length(contentsTdERFStatsFolder)
    if strfind(contentsTdERFStatsFolder(kk).name,'.mat')
        nNeurons = nNeurons + 1;
        neuronNames(nNeurons).name = contentsTdERFStatsFolder(kk).name;
    end
end

%% For each of those neurons, compute the RF, fit a gaussian, plot and save.

fh = figure();

for kk=1:nNeurons
    % Load the ERF movie data
    load(fullfile(tdERFStatsFolder, neuronNames(kk).name));
    cNeuronID = str2double(neuronNames(kk).name(7:end-8));

    plotTimeDependentReceptiveField(neuronERFdata, fh, cNeuronID, outputFolderFigures, nFrames);
end

end % computeElectricReceptiveField

function plotTimeDependentReceptiveField(neuronTdERFdata, fh, neuronID, ...
    outputFolder, maxNumberOfFrames)
% Function dedicated to plotting the eRF for a single neuron.

fileName = fullfile(outputFolder, sprintf('neuron%d_td_erf.gif',neuronID));
nFrames = size(neuronTdERFdata.allPSTHs,2);

% Formatting...
neuronPos = neuronTdERFdata.neuronPos;
allPixelPositions = [neuronTdERFdata.pixelPos; neuronTdERFdata.surroundPixels];
allOutputs = [neuronTdERFdata.allPSTHs; ...
                zeros(size(neuronTdERFdata.surroundPixels,1),nFrames)];
if nnz(allOutputs(:))
    allOutputs = allOutputs/max(allOutputs(:));
end
nPixels = size(neuronTdERFdata.pixelPos,1);

[V,C] = voronoin([allPixelPositions(:,1),allPixelPositions(:,2)]);

% Plotting the result
fh = figure(fh); clf; set(fh,'color','white');

for kk=1:min(nFrames,maxNumberOfFrames)
    clf; hold on;
    
    % Plotting the pixels
    for ii = 1:nPixels
        if all(C{ii}~=1)   % If at least one of the indices is 1, 
                           % then it is an open region and we can't 
                           % patch that.
            patch(V(C{ii},1),V(C{ii},2),repmat(allOutputs(ii,kk),1,3)); % use color i.
        end
    end
    
    % If the neuron position was specified, add a marker
    if ~isempty(neuronPos)
        scatter(neuronPos(1), neuronPos(2), 'b', 'MarkerFaceColor',[0 0 1]);
    end
    
    % Title and stuff
    axis equal; axis off;
    title(sprintf('Neuron %d ERF movie; max #spikes/bin: %2.2f', ...
        neuronID, max(neuronTdERFdata.allPSTHs(:))),...
        'fontsize',14)
    
    % Preparing for saving as gif
    drawnow; 
    frame = getframe(fh);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);

    % Saving
    if kk==1
        imwrite(imind, cm, fileName, 'gif', 'Loopcount', inf, 'DelayTime', 0.2);
    else
        imwrite(imind, cm, fileName, 'gif', 'WriteMode', 'append', 'DelayTime', 0.2);
    end
end

end % plotTimeDependentReceptiveField
