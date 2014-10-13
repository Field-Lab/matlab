function [p, f] = plotElicitedSpikesVsDistance(dataFolder,statsFolder,figuresFolder,spotCenterPosition,varargin)
% plotElicitedSpikesVsDistance(dataFolder,statsFolder,figuresFolder,spotCenterPosition,...)
%
% This function plots the number of spikes elicited by a spot of visible or
% NIR light versus the distance to the center of the spot. It will attempt
% to fit a Gaussian with tails to this plot (note: a supergaussian function
% might be a better model for large spots!), which can be used to estimate
% the spread of the stimulation. 
%
% Parameters:
%   - dataFolder: path to the data folder in which a neurons, params and
%   ei files can be found.
%   - statsFolder: path to the neuron statistics folder.
%   - figuresFolder: path to the folder where the figures will be saved.
%   - spotCenterPosition: position in micrometers of the center of the
%   spot. The point of coordinates (0,0) corresponds to the center of the
%   microelectrode array. See the electrode map for more details on the
%   coordinates.
%
% Additional optional 'parameter'/value pairs:
%   - 'savePlot': by default true, controls if the plots are saved or not.
%   - 'titleString': title of the plots
%   - 'stimulatedNeuronsList': if provided in the form of an array of
%   neuronIDs, the function will only plot data for the neurons in this
%   list.
%   - 'figNumber': figure handle in which the plots will be created. By
%   default fig 1. 
%   - 'stimToPlot': stimulation level that is to be plotted in the dataset.
%   If a dataset consists of k stimuli, the stimToPlot stimulus will be
%   plotted. This index is 1-based, and defaults to 1 (first stimulus).
%
% Returns:
%   - p: parameters of the function that was fit onto the datapoints.
%   - f: the function that was fit onto the datapoints. Typically, 
%   f = @(p,x) p(1).*exp(-abs(x.^p(3))./p(2)^2);
%
% Version: 0.2 - 02/04/2013

%% Reading the input arguments

savePlot = true;
titleString = '#spikes/trial vs. distance to the spot';
figNumber = 1;
powLevelToPlot = 1; % index of the stimulation run
imageFormat = 'epsc';
useSignificanceWeighting = true;
significanceWeightingFunction = @(x) x;
stimulatedNeurons = []; % Default
                 
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
        case 'titlestring'
            titleString = varargin{kk*2};
        case 'fignumber'
            figNumber = varargin{kk*2};
        case 'stimtoplot'
            powLevelToPlot = varargin{kk*2};
        case 'imageformat'
            imageFormat = varargin{kk*2};
        case 'stimulatedneuronslist'
            stimulatedNeurons = varargin{kk*2};
        case 'saveplot'
            savePlot = varargin{kk*2};
        otherwise
            err = MException('MATLAB:InvArgIn',...
                'Unknown parameter specified');
            throw(err);
    end
end

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
if size(spotCenterPosition,2) > size(spotCenterPosition,1)
    spotCenterPosition = spotCenterPosition.';
end

if ~exist(figuresFolder,'dir')
    mkdir(figuresFolder);
end


%% Some default values set here

if isunix
    visionPath = '/home/ggoetz/Research/Vision/Vision815/Vision.jar';
else
    visionPath = 'C:\Users\ggoetz\Research\Vision\Vision815\Vision.jar';
end

%% Linking to the parameter and neuron files

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

%% From param file, compute distance between each neuron EI and spot center

idList=neuronFile.getIDList();
if ~isempty(stimulatedNeurons)
    idList = intersect(idList, stimulatedNeurons);
end
if size(idList,1)<size(idList,2)
    idList = idList.';
end
nNeurons = length(idList);
distToSpot = zeros(nNeurons,1);

% Computing the distances to the spot
for ii=1:nNeurons
    eiPos = [paramFile.getDoubleCell(idList(ii), 'EIx0'); ...
             paramFile.getDoubleCell(idList(ii), 'EIy0')];
    distToSpot(ii) = norm(spotCenterPosition-eiPos);
end

%% Read the activation data

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

% Loading the activation data in the second column of idList
actData = zeros(size(idList,1),2)*NaN;

load([statsFolder neuronNames(1).name]);
stimDataPos = ismember(activationData.labels,'nSpikes');
errorDataPos = ismember(activationData.labels,'significance');
if sum(errorDataPos)>0
    useSignificanceWeighting = true && useSignificanceWeighting;
else
    useSignificanceWeighting = false;
end

for kk=1:nNeurons
    pos = idList == neuronNames(kk).ID;
    load([statsFolder neuronNames(kk).name]);
    actData(pos,1) = activationData.data(powLevelToPlot,stimDataPos);
    if useSignificanceWeighting
        actData(pos,2) = activationData.data(powLevelToPlot,errorDataPos);
    else
        actData(pos,2) = 1;
    end
end

% Removing the neurons for which the distance data or the activation data
% is incorrect

pos_incorrect = logical(isnan(actData(:,1)) + isnan(actData(:,2)) + isnan(distToSpot));
actData(pos_incorrect,:) = [];
distToSpot(pos_incorrect) = [];

%% Plot the result

% Plotting the neuron response vs distance to the spot
fh = figure(figNumber); clf; 
set(fh,'color','white'); hold on;
scatter(distToSpot,actData(:,1),'.k');

% Fit a (quasi) Gaussian onto the datapoints
% f = @(p,x) p(1).*exp(-x.^2./(2*p(2)^2));
f = @(p,x) p(1).*exp(-abs(x.^p(3))./p(2)^2);

% Specifying fitting options
fo_ = fitoptions('method','NonlinearLeastSquares',...
    'Lower',[0 10 1.5],'Upper',[max(1.5*max(actData(:,1)),0.1) 1000 2]);
st_ = [10 200 1.5];
set(fo_,'Startpoint',st_);
if useSignificanceWeighting
    weightData = significanceWeightingFunction(actData(:,2));
    weightData(weightData==Inf) = max((weightData(weightData~=Inf))*10^2); 
    set(fo_,'Weight',weightData);
end
ft_ = fittype('a*exp(-abs(x^c)/b^2)',...
    'dependent',{'y'},'independent',{'x'},...
    'coefficients',{'a', 'b', 'c'});

% Doing the fitting
cf_ = fit(distToSpot,actData(:,1),ft_,fo_);
p = coeffvalues(cf_);

% Add to plot
line(linspace(0,max(distToSpot),400),f(p,linspace(0,max(distToSpot),400)),...
    'color','r')

% Legend, title and stuff
title(titleString,'fontsize',14)
grid on
xlabel('Distance to the spot, \mum');
ylabel('#spikes/trial');
axis([0 1000 1.1*min(0,min(actData(:,1))) 1.1*max(max(actData(:,1)),f(p,0))]);

%% Saving

if savePlot
    saveas(fh,[figuresFolder 'stimVsDist' num2str(powLevelToPlot-1)],imageFormat);
end

end % plotElicitedSpikesVsDistance