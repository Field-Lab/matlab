function plotPositionEI(dataFolder,statsFolder,figuresFolder,varargin)
% plotPositionEI(dataFolder,statsFolder,figuresFolder,...)
% 
% Parameters:
%   - dataFolder:
%   - statsFolder:
%   - figuresFolder:
%
% Optional parameters/value pairs:
%   - displayNeuronNames: by default false. If set to true, the neuron
%   numbers will be plotted on the map.
%   - titleString: by default 'Activation map', title of the figure. 
%   - figNumber: by default 1. Handle to the figure in which the map will
%   be plotted.
%   - minNSpikes: maximum number of spikes per trial over which the color
%   does not change anymore on the map. By default 2.5.
%   - maxNSpikes: minimum number of spikes per trial under which the color
%   doesn't change anymore on the map. 
%   - stimToPlot: by default the functions will plot the activation map for
%   the first stimulus among those ran for this experiment. Using this
%   parameter the user can specify which stimulus should be plotted.
%   - imageFormat: by default .png. Other supported formats include (and
%   are not limited to) .eps, .jpg, .bmp.
%   - radiusPixel: by default 70um, corresponds to the diameter of the
%   circle plotted, whose center is specified in positionPixel.
%   - positionPixel: by default empty. If specified, it should be given in
%   the form of the electrode(s) number(s) closest to the center(s) of the
%   pixel(s)/spot(s) of light.
%   - stimulatedNeuronsList: by default this list is empty, in which case
%   the activation for all neurons will be plotted. If this list is
%   specified, only the neurons in the list will be plotted. Neurons
%   specified but not found in the neurons file are ignored.
%   - savePlot: default is true. Will save the plot in the figuresFolder
%   folder using the default or specified image format.
%   - nColorsColormap: number of different colors in the colormap.
%
% Returns:
%   []
%
% Version: 0.3 - 07/06/2012
%

displayNeuronNames = false;
savePlot = true;
titleString = 'Activation map';
figNumber = 1;
minNSpikes = 0;
maxNSpikes = 2.5;
powLevelToPlot = 1; % index of the stimulation run
imageFormat = 'png';
nColorsColormap = 10;

posPixels = []; % Electrode under the center of each pixel
radPixels = 70;

stimulatedNeurons = []; % Default
% stimulatedNeurons = [63 242 376 499 528 663 812 902 1007 1037 1096 1263 ...
%                      1306 1352 1486 1563 1577 1652 1877 2088 2117 2147 ...
%                      2206 2343 2492 2671 2731 3181 3211 3557 3587 ...
%                      3872 4157 4396 4427 4636 4653 4741 ...
%                      5012 5086 5087 5116 5177 5311 5387 5538 ...
%                      5523 5671 5747 6033 6046 6047 6512 6873 6917 ...
%                      7082 7231].'; % 2012-04-19-0 list

                 
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
        case 'displayneuronnames'
            displayNeuronNames = varargin{kk*2};
        case 'titlestring'
            titleString = varargin{kk*2};
        case 'fignumber'
            figNumber = varargin{kk*2};
        case 'minnspikes'
            minNSpikes = varargin{kk*2};
        case 'maxnspikes'
            maxNSpikes = varargin{kk*2};
        case 'stimtoplot'
            powLevelToPlot = varargin{kk*2};
        case 'imageformat'
            imageFormat = varargin{kk*2};
        case 'radiuspixel'
            radPixels = varargin{kk*2};
        case 'positionpixel'
            posPixels = varargin{kk*2};
        case 'stimulatedneuronslist'
            stimulatedNeurons = varargin{kk*2};
        case 'saveplot'
            savePlot = varargin{kk*2};
        case 'ncolorscolormap'
            nColorsColormap = varargin{kk*2};
        otherwise
            err = MException('MATLAB:InvArgIn',...
                'Unknown parameter specified');
            throw(err);
    end
end

%% Default values
if isunix
    javaaddpath '/home/ggoetz/Research/Vision/Vision815/Vision.jar';
    % Default path to the electrodes coordinates file
    electrodeMapPath = './512coords.txt';
else
    javaaddpath 'C:\Users\ggoetz\Research\Vision\Vision815\Vision.jar';
    % Default path to the electrodes coordinates file
    electrodeMapPath = '.\512coords.txt';
end

% Offsets for positioning the electrode map
mapXOffset = 1042;
mapYOffset = 845;

% Default plot colors for different neuron groups
displayType.marker = ['o';'.';'+';'^';'*';'x';'^';'+']; 
displayType.color = ['b';'r';'g';'c';'m';'b';'r';'g'];


%% Reading the input arguments

% Making sure folders ends by '\' or '/'
if dataFolder(end:end)~=filesep
    dataFolder = [dataFolder filesep];
end
if statsFolder(end:end)~=filesep
    statsFolder = [statsFolder filesep];
end
if figuresFolder(end:end)~=filesep
    figuresFolder = [figuresFolder filesep];
end

if ~exist(figuresFolder,'dir')
    mkdir(figuresFolder);
end

%% Linking to the files

if ~exist('edu/ucsc/neurobiology/vision/io/NeuronFile','class')
    javaaddpath(visionPath);
end

% Finding the neuron file, the param file and the ei file
contentsDataFolder = dir(dataFolder);
for kk=1:length(contentsDataFolder)
    isNeuronFile = strfind(contentsDataFolder(kk).name,'.neurons');
    isNeuronRawFile = strfind(contentsDataFolder(kk).name,'.neurons-raw');
    isParamFile = strfind(contentsDataFolder(kk).name,'.params');
    if isNeuronFile
        if isNeuronRawFile
        else
            if (isNeuronFile+length('.neurons')-1)==length(contentsDataFolder(kk).name)
                neuron_path =  [dataFolder contentsDataFolder(kk).name];
            end
        end
    end
    if isParamFile
        param_path = [dataFolder contentsDataFolder(kk).name];
    end
end
% Linking to the files
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile(neuron_path);
paramFile=edu.ucsc.neurobiology.vision.io.ParametersFile(param_path);

%% Reading the params files to get EI position

idList=neuronFile.getIDList();
if ~isempty(stimulatedNeurons)
    idList = intersect(idList, stimulatedNeurons);
end
if size(idList,1)<size(idList,2)
    idList = idList.';
end
nNeurons = length(idList);
eiPosAll = zeros(nNeurons,2);

% Reading the EI positions - warning: some of them could be NaN
for ii=1:nNeurons
    eiPosAll(ii,1)=paramFile.getDoubleCell(idList(ii), 'EIx0');
    eiPosAll(ii,2)=paramFile.getDoubleCell(idList(ii), 'EIy0');
end

%% Reading the electrode map

% Note: 512coords.txt has the coordinates of the electrodes in microns
posElectrodes = dlmread(electrodeMapPath);   

%% Reading the activation data 

% Finding the activation files
contentsStatsFolder = dir(statsFolder);
neuronNames = struct('name','');
nNeurons = 0;

for kk=1:length(contentsStatsFolder)
    pos = strfind(contentsStatsFolder(kk).name,'.mat');
    if (pos)
        nNeurons = nNeurons + 1;
        neuronNames(nNeurons).name = contentsStatsFolder(kk).name;
        neuronNames(nNeurons).ID = str2double(contentsStatsFolder(kk).name(7:(pos-5)));
    end
end

% Loading the activation data
idList = double([idList zeros(size(idList))]);

load([statsFolder neuronNames(1).name]);
stimDataPos = ismember(activationData.labels,'nSpikes');
for kk=1:nNeurons
    pos = idList(:,1) == neuronNames(kk).ID;
    load([statsFolder neuronNames(kk).name]);
%     idList(pos,2) = max(activationData.data(:,stimDataPos));
    idList(pos,2) = activationData.data(powLevelToPlot,stimDataPos);
end

nNeurons = length(idList);

%% Plotting

% Plotting the electrode map
fh = figure(figNumber); clf; 
set(fh,'color','white');
set(fh, 'Position', [100 100 1280 720])
scatter((posElectrodes(:,1) + mapXOffset)*2,(posElectrodes(:,2) + mapYOffset)*2,'.k');
hold on

% Plotting the pixels
for kk=1:length(posPixels)
    circle([(posElectrodes(posPixels(kk)+1,1) + mapXOffset)*2,...
        (posElectrodes(posPixels(kk)+1,2) + mapYOffset)*2],...
        radPixels*2,1000,'r');
end

axis equal;

% Getting the colormap for the neurons
colormap(flipud(hot(nColorsColormap)));

% Plotting the neurons
for kk=1:nNeurons
    cActivation = idList(kk,2);
    
    scatter((eiPosAll(kk,1) + mapXOffset)*2, ...
            (eiPosAll(kk,2) + mapYOffset)*2,...
            60 + 20*max(cActivation,0),cActivation,'o','filled');

    if displayNeuronNames
        text((eiPosAll(kk,1) + mapXOffset)*2, ...
                (eiPosAll(kk,2) + mapYOffset)*2,...
                ['  ',num2str(idList(kk))], 'FontSize',8);
    end

end

hold off

% Colorbar, legend and stuff
axis off
title(titleString,'fontsize',14)
t = colorbar();
set(get(t,'ylabel'),'string','#spikes elicited','Fontsize',10)
set(t,'position',[0.92 0.22 0.023 0.65])
caxis([minNSpikes maxNSpikes])

%% Saving

if savePlot
    saveas(fh,[figuresFolder 'activation_stim' num2str(powLevelToPlot-1)],imageFormat);
end

end

function H=circle(center,radius,NOP,style)
%---------------------------------------------------------------------------------------------
% H=CIRCLE(CENTER,RADIUS,NOP,STYLE)
% This routine draws a circle with center defined as
% a vector CENTER, radius as a scaler RADIS. NOP is 
% the number of points on the circle. As to STYLE,
% use it the same way as you use the rountine PLOT.
% Since the handle of the object is returned, you
% use routine SET to get the best result.
%
%   Usage Examples,
%
%   circle([1,3],3,1000,':'); 
%   circle([2,4],2,1000,'--');
%
%   Zhenhai Wang <zhenhai@ieee.org>
%   Version 1.00
%   December, 2002
%---------------------------------------------------------------------------------------------

if (nargin <3),
 error('Please see help for INPUT DATA.');
elseif (nargin==3)
    style='b-';
end;
THETA=linspace(0,2*pi,NOP);
RHO=ones(1,NOP)*radius;
[X,Y] = pol2cart(THETA,RHO);
X=X+center(1);
Y=Y+center(2);
H=plot(X,Y,style);
axis square;

end % circle
