function allThresholds = plotActivation(statsFolder,figuresFolder,thresholdsFolder,varargin)
% TODO: documentation
% 
% Function used to plot sigmoid activation curves. 
%
% Parameters:
%   - statsFolder: path to the folder containing the activation statistics
%   (not the stats data!).
%   - figuresFolder: path to the folder containing the activation plots
%   - thresholdsFolder: path to the folder where the thresholds should be
%   saved.
%
% Optional value/parameters combinations:
%   - plotType: see description of the variable in the code. For irradiance
%   thresholds, use 1 (default value).
%   - fitCurve
%   - findThreshold
%   - imageFormat
%
% Returns:
%    - allThresholds: matrix of all the thresholds found. First column is
%    the neuronID, second column the chi^2 of the last fit and the
%    following columns the thresholds.
%
% Version: v4.06 - 05/29/2013
%

% TODO: change chi^2 storage: there should be one chi^2 per fit

%% Reading the input arguments, setting default values

fh = 0;

plotType = 1;   % 1 = different irradiances on same plot; fixed PW, frequency
                % 2 = different PW on same plot; fixed irradiance, frequency
                % 3 = different frequencies on same plot; fixed PW, irradiance
                % 4 = different irradiances on same plot, fixed image, frequency
                
pauseBetweenPlots = false;
fitCurve = true;
findThreshold = true;
savePlots = true;
% displayPlots = true;        % Unused for now
saveThresholds = true;

maxFigsAuthorized = 1;
imageFormat = 'png';

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
        case 'plottype'
            plotType = varargin{kk*2};
        case 'fitcurve'
            fitCurve = varargin{kk*2};
        case 'findthreshold'
            findThreshold = varargin{kk*2};
        case 'imageformat'
            imageFormat = varargin{kk*2};
        otherwise
            err = MException('MATLAB:InvArgIn',...
                'Unknown parameter specified');
            throw(err);
    end
end

% Formatting the folders
if statsFolder(end:end)~=filesep
    statsFolder = [statsFolder filesep];
end
if figuresFolder(end:end)~=filesep
    figuresFolder = [figuresFolder filesep];
end
if thresholdsFolder(end:end)~=filesep
    thresholdsFolder = [thresholdsFolder filesep];
end

%% Checking folders, getting the name of all the neurons

if ~exist(figuresFolder,'dir')
    mkdir(figuresFolder);
end
if ~exist(thresholdsFolder,'dir')
    mkdir(thresholdsFolder);
end

contentsStatsFolder = dir(statsFolder);
neuronNames = struct('name','');
nNeurons = 0;

for kk=1:length(contentsStatsFolder)
    if strfind(contentsStatsFolder(kk).name,'.mat')
        nNeurons = nNeurons + 1;
        neuronNames(nNeurons).name = contentsStatsFolder(kk).name;
    end
end

%% Figuring out how many plots we're going to have
% Create the thresholds matrix accordingly

load([statsFolder neuronNames(1).name]);
[~, nCombinations] = formatData(plotType,activationData);

% First column: neuronID
allThresholds = NaN*zeros(nNeurons,nCombinations+1);


%% Getting the plot with error bars
for kk=1:nNeurons
    %% Getting the experimental data
    load([statsFolder neuronNames(kk).name]);
    display(['Plotting data for neuron ' neuronNames(kk).name(7:end-8)])
    allThresholds(kk,1) = str2double(neuronNames(kk).name(7:end-8));
    
    % Defining the points of interest
    [plotData, nCombinations] = formatData(plotType,activationData);
    sizeSubFig = ceil(sqrt(nCombinations));
    
    %% Creating the plots
    if nNeurons<=maxFigsAuthorized
        fh = figure(fh+1); clf; set(fh,'color','white');
    else
        fh = figure(max(fh,1)); clf; set(fh,'color','white');
    end
    set(fh,'position',[100 100 800 600])

    for ll=1:nCombinations
        % Selecting the subfigure in which the plot goes
        subplot(sizeSubFig,sizeSubFig,ll)
        
        % Selecting the data from the plotData structure
        y = plotData(ll).stimData; 
        x = plotData(ll).varyingParam;
        e = plotData(ll).errorData;

        errorbar(x,y,sqrt(e),'MarkerFaceColor',[0 0 0],'Marker','square',...
            'LineStyle','none',...
            'Color',[0 0 0]);
        switch plotType
            case 1
                xlabel('Irradiance (mW/mm^2)')
                legendStr = [num2str(plotData(ll).fixedParam,...
                    '%1.1f') plotData(ll).fixedParamUnit{1}];
                title(['Pulse duration: ' legendStr]);
            case 2
                xlabel('Pulse duration (ms)')
                legendStr = [num2str(plotData(ll).fixedParam,...
                    '%1.3f') plotData(ll).fixedParamUnit{1}];
                title(['Irradiance: ' legendStr])
            case 3
                xlabel('Pulse frequency (Hz)')
                legendStr = [num2str(plotData(ll).fixedParam,...
                    '%1.1f') plotData(ll).fixedParamUnit{1}];
                title(['Pulse duration: ' legendStr]);
            case 4
                xlabel('Irradiance (mW/mm^2)')
                legendStr = ['Image ' num2str(plotData(ll).fixedParam)];
                title(legendStr);
        end
        ylabel('#spikes/trial')
        axis([0 max(x)*1.05 min(0,min(y)*1.05) max(0.5,max(y)*1.05)])
        
        if fitCurve
            % Trying to fit the cdf of a gamma distribution
            % See curve fitting toolbox for explanation of the code
            f = @(p,x) p(1)*gamcdf(x,p(2),p(3));
            w = 1./e;
            % Infinite weights should be removed
            if sum(w~=Inf)
                w(w==Inf) = max((w(w~=Inf))*2); 
            else
                w(w==Inf) = 1;
            end
            
            % Specifying fitting options
            fo_ = fitoptions('method','NonlinearLeastSquares',...
                'Lower',[0 0.05 0.05],'Upper',[10 10 10]);
            st_ = [1 0.5 1];
            set(fo_,'Startpoint',st_);
            set(fo_,'Weight',w);
            
            ft_ = fittype('a*gamcdf(x,b,c)',...
                'dependent',{'y'},'independent',{'x'},...
                'coefficients',{'a', 'b', 'c'});

            % Doing the fitting
            cf_ = fit(x,y,ft_,fo_);
            p = coeffvalues(cf_);
            
%             % Computing the reduced chi-square
%             allThresholds(kk,2) = sum((y-f(p,x)).^2./e);

            %  Adding the fit to the plot
            line(linspace(-0,max(x),400),f(p,linspace(-0,max(x),400)),...
                'color','r')
            
            % If we did fitting we can go look for a threshold
            if findThreshold
                % Computing the threshold
                f_thresh = @(x) p(1)*gamcdf(x,p(2),p(3))-0.5;
                thresh = fzero(f_thresh,mean(x),optimset('Display','off'));
                
                % Discarding thresholds found outside of the measured range
                if (thresh>max(x))||(thresh<min(x))
                    thresh = NaN;
                end
                
                % Displaying and storing it
                if ~isnan(thresh)
                    allThresholds(kk,ll+1) = thresh;
                    display(['   ' legendStr ' pulse: threshold = '...
                        num2str(thresh,'%2.4f')])
                else
                    display(['   ' legendStr ' pulse: no threshold found'])
                end
            end
        end
            
    end
    
    % Adding a title at the top
    axes('Position', [0 0 1 1], 'Visible', 'Off');
    h = text(0.5,0.99, ['Neuron ' neuronNames(kk).name(7:end-8)],...
        'HorizontalAlignment','Center','VerticalAlignment','Top');
    set(h,'FontSize',16)
    
    if savePlots
        saveas(fh,[figuresFolder neuronNames(kk).name(1:end-8)],imageFormat);
    end
    
    if pauseBetweenPlots
        pause
    end
end

if saveThresholds
    save([thresholdsFolder 'allThresholds.mat'],'allThresholds');
end

end % plotStrDur

function [plotData, nCombinations] = formatData(plotType, activationDataStruct)
% This function formats the activation data in a way that makes it usable
% for subsequently plotting strength/duration curves.

activationData = activationDataStruct.data;

% Finding where the data is stored
activationLabels = activationDataStruct.labels;
stimDataPos = ismember(activationLabels,'nSpikes');
errorDataPos = ismember(activationLabels,'SEM');
irradianceDataPos = ismember(activationLabels,'Irradiance-mW/mm^2');
PWDataPos = ismember(activationLabels,'Duration');
freqDataPos = ismember(activationLabels,'Frequency');
imDataPos = ismember(activationLabels,'Image projected');

% Building the plot Data structure
plotData = struct('pos',{},'stimData',{},'errorData',{},...
    'varyingParam',{},'fixedParam',{},'fixedParamUnit',{});
units = {'ms','mW/mm^2','Hz'};

% Finding where the fixed parameter column is and where the varying
% parameter column is
switch plotType
    case 1
        fixedParamPos = PWDataPos;
        varyingParamPos = irradianceDataPos;
        unitPos = 1;
    case 2
        fixedParamPos = irradianceDataPos;
        varyingParamPos = PWDataPos;
        unitPos = 2;
    case 3
        fixedParamPos = PWDataPos;
        varyingParamPos = freqDataPos;
        unitPos = 1;
    case 4
        fixedParamPos = imDataPos;
        varyingParamPos = irradianceDataPos;
        unitPos = 1;
    otherwise
        err = MException('MATLAB:InvArgIn',...
            'Unknown plot type specified');
        throw(err);
end

% Getting all the different values for the fixed parameter
allFixedParams = unique(activationData(:,fixedParamPos));
nCombinations = length(allFixedParams);
for kk=1:nCombinations
    plotData(kk).fixedParamUnit = units(unitPos);
end

% For each one of these values we 
for kk=1:nCombinations
    plotData(kk).fixedParam = allFixedParams(kk);
    plotData(kk).pos = find(activationData(:,fixedParamPos)==allFixedParams(kk));
    
    plotData(kk).varyingParam = activationData(plotData(kk).pos,varyingParamPos);
    plotData(kk).stimData = activationData(plotData(kk).pos,stimDataPos);
    plotData(kk).errorData = activationData(plotData(kk).pos,errorDataPos);
end

end % findPOI

function y = genLogFcn(p,x)

% Renaming the variables for coherence with Wikipedia article on the 
% generalized logistic function
A = 0;
K = p(1);
B = p(2);
Q = p(3);
nu = p(4);
M = p(5);

y = A + (K-A)./((1 + Q.*exp(-B.*(x-M))).^(1./nu));

end