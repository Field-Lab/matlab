function plotElicitedSpikesVsDistance2Px(dataFolder,statsFolder,figuresFolder,spotCenterPosition,varargin)
% plotElicitedSpikesVsDistance(dataFolder,statsFolder,figuresFolder,spotCenterPosition,...)
%
% Parameters:
%   - dataFolder:
%   - statsFolder:
%   - figuresFolder:
%   - spotCenterPosition: it should be a 4*1 or 1*4 matrix such that the
%   positions of the first pixels are [x,y] = spotCenterPosition(1:2) and
%   for the second pixel, [x,y] = spotCenterPosition(3:4).
%
% Optional parameters
%   - titleString: title of the plot
%   - figNumber: number of the figure in which things will be plotted.
%   - stimToPlot: stimulation level to plot. 1 is the first stimulus, 2 the
%   second one and so on.
%   - imageFormat: string, specifies in what format the image will be
%   saved.
%   - stimulatedNeuronsList: list of the neurons that were stimulated. If
%   this list is empty, all the neurons in the neurons file will be
%   included in the plot.
%   - savePlot: by default true, if set to false the plot will not be
%   saved.
%
% Returns:
%
% Version: 1.0 - 03/26/2013

%% Reading the input arguments

savePlot = true;
titleString = '';
figNumber = 1;
powLevelToPlot = 1; % index of the stimulation run
imageFormat = 'png';
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
eiPos = zeros(nNeurons,2);
distToSpot = zeros(nNeurons,2);

% Computing the distances to the spot
for ii=1:nNeurons
    eiPos(ii,:) = [paramFile.getDoubleCell(idList(ii), 'EIx0'); ...
             paramFile.getDoubleCell(idList(ii), 'EIy0')];
    distToSpot(ii,1) = norm(spotCenterPosition(1:2).'-eiPos(ii,:));
    distToSpot(ii,2) = norm(spotCenterPosition(3:4).'-eiPos(ii,:));
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

pos_incorrect = logical(isnan(actData(:,1)) + isnan(actData(:,2)) + ...
    isnan(eiPos(:,1)) + isnan(eiPos(:,2)));
actData(pos_incorrect,:) = [];
eiPos(pos_incorrect) = [];
distToSpot(pos_incorrect) = [];

%% Plot the result

%%% Plot 1: MEA fit

fh = figure(figNumber); clf; 
set(fh,'color','white'); hold on;

% Plotting the EI map
scatter(eiPos(:,1),eiPos(:,2),50*ones(size(eiPos(:,2))),actData(:,1),'filled');

% Fit a sum of two 2D exponential functions onto the data
f = @(p,x,y,spotCenterPosition) p(1).*exp(-abs(((x-spotCenterPosition(1)).^p(2)+(y-spotCenterPosition(2)).^p(2)))./p(3)^2) + ...
    p(4).*exp(-abs(((x-spotCenterPosition(3)).^p(5)+(y-spotCenterPosition(4)).^p(5)))./p(6)^2);

% Specifying fitting options
fo_ = fitoptions('method','NonlinearLeastSquares',...
    'Lower',[0 1.5 10 0 1.5 10],'Upper',[15 2 1000 15 2 1000]);
st_ = [10 1.5 100 10 1.5 100];
set(fo_,'Startpoint',st_);
if useSignificanceWeighting
    weightData = significanceWeightingFunction(actData(:,2));
    weightData(weightData==Inf) = max((weightData(weightData~=Inf))*10^2); 
    set(fo_,'Weight',weightData);
end
ft_ = fittype(@(a,b,c,d,e,f,x,y) a.*exp(-abs(((x-spotCenterPosition(1)).^b+(y-spotCenterPosition(2)).^b))./c^2) + ...
    d.*exp(-abs(((x-spotCenterPosition(3)).^e+(y-spotCenterPosition(4)).^e))./f^2),...
    'dependent',{'z'},'independent',{'x','y'},...
    'coefficients',{'a', 'b', 'c', 'd', 'e', 'f'});

% Doing the fitting
cf_ = fit(eiPos,actData(:,1),ft_,fo_);
p = coeffvalues(cf_);

% Add to plot
f_ez = @(x,y) f(p,x,y,spotCenterPosition);
ezcontour(f_ez,[1.1*min(eiPos(:,1)) 1.1*max(eiPos(:,1)) ...
    1.1*min(eiPos(:,1)) 1.1*max(eiPos(:,1))],150);

% Legend, title and stuff
title(titleString,'fontsize',14)
grid on
xlabel('x, \mum');
ylabel('y, \mum');
axis([-900 900 -500 500]); axis equal; 
t = colorbar();
set(get(t,'ylabel'),'string','#spikes elicited','Fontsize',10)

%%% Plot 2: 3D surface plot

fh = figure(fh+1); clf; 
set(fh,'color','white'); %,'renderer','openGL');

ezsurf(f_ez,[1.1*min(eiPos(:,1)) 1.1*max(eiPos(:,1)) ...
    1.1*min(eiPos(:,1)) 1.1*max(eiPos(:,1))],40);
alpha(0.9);
hold on;
scatter3(eiPos(:,1),eiPos(:,2),actData(:,1),'k','filled');
% scatter3(eiPos(:,1),eiPos(:,2),actData(:,1),50*ones(size(eiPos(:,2))),actData(:,1),'filled');

% axis([-900 900 -900 900 0 1.1*max(actData(:,1))])

title(titleString,'fontsize',14)
grid on
xlabel('x, \mum');
ylabel('y, \mum');
zlabel('#spikes/trial');
% axis equal; % axis([-900 900 -500 500]); 
t = colorbar();
set(get(t,'ylabel'),'string','#spikes elicited','Fontsize',10)

%%% Plot 3: deviation from fit

fh = figure(fh+1); clf; 
set(fh,'color','white'); hold on;

deviationFromFit = actData(:,1) - f_ez(eiPos(:,1),eiPos(:,2));
scatter(eiPos(:,1),eiPos(:,2),50*ones(size(eiPos(:,2))),deviationFromFit(:,1),'filled');
axis([-900 900 -500 500]); axis equal; 
t = colorbar();
set(get(t,'ylabel'),'string','deviation from fit in #spikes','Fontsize',10)


%%% Plot 4: #spikes vs dist to pixels

distBetweenPixels = norm(spotCenterPosition(1:2)-spotCenterPosition(3:4));

% Fit parameters here
f = @(p,x,y,distBetweenPixels) p(1).*exp(-abs(((x-distBetweenPixels).^p(2)+y.^p(2)))./p(3)^2) + ...
    p(4).*exp(-abs((x.^p(5)+(y-distBetweenPixels).^p(5)))./p(6)^2);
fo_ = fitoptions('method','NonlinearLeastSquares',...
    'Lower',[0 1.5 10 0 1.5 10],'Upper',[15 2 1000 15 2 1000]);
st_ = [10 1.5 100 10 1.5 100];
set(fo_,'Startpoint',st_);
if useSignificanceWeighting
    weightData = significanceWeightingFunction(actData(:,2));
    weightData(weightData==Inf) = max((weightData(weightData~=Inf))*10^2); 
    set(fo_,'Weight',weightData);
end
ft_ = fittype(@(a,b,c,d,e,f,x,y) a.*exp(-abs(((x-distBetweenPixels).^b+y.^b))./c^2) + ...
    d.*exp(-abs((x.^e+(y-distBetweenPixels).^e))./f^2),...
    'dependent',{'z'},'independent',{'x','y'},...
    'coefficients',{'a', 'b', 'c', 'd', 'e', 'f'});

% Doing the fitting
cf_ = fit(distToSpot,actData(:,1),ft_,fo_);
p = coeffvalues(cf_);

% Plot
fh = figure(fh+1); clf; 
set(fh,'color','white'); hold on;
scatter(distToSpot(:,1),distToSpot(:,2),50*ones(size(eiPos(:,2))),actData(:,1),'filled');
f_ez = @(x,y) f(p,x,y,distBetweenPixels);
ezsurfc(f_ez,[0 1 0 1]*max(max(distToSpot))*1.1,50)

% Title and stuff
title(titleString,'fontsize',14)
grid on
xlabel('Distance to spot 1, \mum');
ylabel('Distance to spot 2, \mum');
% axis([-900 900 -500 500]); axis equal; 
t = colorbar();
set(get(t,'ylabel'),'string','#spikes elicited','Fontsize',10)

%% Saving

fh = figNumber;
if savePlot
    saveas(fh,[figuresFolder 'stimVsDist2D_plot1_' num2str(powLevelToPlot-1)],imageFormat);
    saveas(fh+1,[figuresFolder 'stimVsDist2D_plot2_' num2str(powLevelToPlot-1)],imageFormat);
    saveas(fh+2,[figuresFolder 'stimVsDist2D_plot3_' num2str(powLevelToPlot-1)],imageFormat);
    saveas(fh+3,[figuresFolder 'stimVsDist2D_plot4_' num2str(powLevelToPlot-1)],imageFormat);
end

end % plotElicitedSpikesVsDistance