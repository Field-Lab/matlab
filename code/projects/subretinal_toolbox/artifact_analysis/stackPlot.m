%%% Stack Plot
% Version: v1.20 - 02/08/2014 - Philip Haeusser
% 

function varargout = stackPlot(varargin)

clear all;

%% Preliminaries

% Specify path to stored data files - include terminal backslash everywhere
dataRoot = '/Users/Philip/UCSC/Retina/artifactAnalysis/';
dataDir = 'run3/';  % e.g. 'data/'
dataPath = strcat(dataRoot,dataDir);
processedDataPath = '/Users/Philip/UCSC/Retina/data/2014-01-17-0/data/processed_run3/data000/';
avgArtifactFile = 'avArtifact006.mat'; % e.g. 'avArtifact006.mat'

% Specify how many electrode files are plotted, starting with
% lowest electrodeRef (not including electrod 0 of course)
howManyElectrodes = 1; % (0 = all)

% Specify how many pulses are plotted
howManyPulses = 0; % (0 = all)

% Duration of pulse in samples (only this many samples will be plotted)
pulseDuration = 1000;

% Left offset (plot this additional part before pulse)
lOffset = 79;

% Vertical offset
vOffset = 300;

% Skip samples
skipSamples = 0;

% Y axis max/min (comment out to set to auto)
% yAxMin = -2500;
% yAxMax = 13500;

% Plot colors?
plotColors = false;

% Plot average artifact??
plotAvArtifact = true;

% Artifact y scaling factor
aScaleFactor=1/2;

% Processed data y scaling factor
pScaleFactor=0.01; %3

% close all figures after plotting?
closeAllFigures = 0;

% save plots to file?
savePlot = 0;

%% Prepare Matlab Garbage Collection
% To call, execute: org.dt.matlab.utilities.JavaMemoryCleaner.clear(0);

MatlabGCPath = which('MatlabGarbageCollector.jar');
if ~exist('org.dt.matlab.utilities.JavaMemoryCleaner','class')
    javaaddpath(MatlabGCPath);
end

%% Get File list
matFiles = dir(strcat(dataPath,'*.mat'));
electrodeZero = matFiles(1).name;
    if ~electrodeZero == 'electrode_0.mat'
        error('Data Folder doesn´t contain electrode_0.mat!');
    end

%% Find pulse times
load(strcat(dataPath,electrodeZero));
pulseDataGrad = abs(diff(originalData));
pulseTimes = find(pulseDataGrad>=.25*max(pulseDataGrad));

% Note: In this pulseTimes list, every other entry should be discarded.
                
%% Load and plot data sets

% Check if all electrode files should be plotted or not
maxElectrode=length(matFiles);
if ~howManyElectrodes == 0;
    maxElectrode=howManyElectrodes+1;
end

% Check if all pulses should be plotted or not
maxPulse=size(pulseTimes,1)/2-3;
if ~howManyPulses == 0;
    maxPulse=min(howManyPulses,size(pulseTimes,1)/2-3);
end

for el=2:maxElectrode;
    % Load original and processed data
    load(strcat(dataPath,matFiles(el).name));
    
    % Get some fancy colors (if desired)
    if plotColors;
        colors=hsv(size(pulseTimes,1)/2);
    end
    
    % Plot intervals
    myPlot = figure; 
    hold on;
    
    % Don't look too far to the left/right (avoid Index exceeds matrix
    % dimensions error)
    lOffset = min(pulseTimes(1)-1, lOffset);
    
    processedData(:)=processedData(:).*pScaleFactor;
    
    vOffset = mean(processedData(...
                pulseTimes(2*skipSamples+1)-lOffset:...
                min(pulseTimes(2*skipSamples+1)+pulseDuration,size(processedData,1))));
    
    for i=0+skipSamples:maxPulse;
        if plotColors
            plot(processedData(...
                max(pulseTimes(2*i+1)-lOffset,1):min(pulseTimes(2*i+1)+pulseDuration,size(processedData)))...
                +vOffset*i,...
                'color',colors(i+1,:));
        else
            plot(processedData(...
                pulseTimes(2*i+1)-lOffset:...
                min(pulseTimes(2*i+1)+pulseDuration,size(processedData,1)))...%.*pScaleFactor...
                -vOffset...
                +i+1);
        end
    end
    
    if plotAvArtifact;
        % align artifact
        load(strcat(processedDataPath,'artifact/',avgArtifactFile));
        
        % get "time zero" 
        avgADataGrad = abs(diff(avArtifact(:,electrodeRef+1)));
        artifactTimes = find(avgADataGrad>=.25*max(avgADataGrad));
        artifactTimeZero = artifactTimes(1);
        avAfShift = artifactTimeZero-lOffset;
        
        % vertical shift
%         zeroLine = mean(processedData(...
%                 max(pulseTimes(1)-lOffset,1):...
%                 min(pulseTimes(1)+pulseDuration,size(processedData,1))))...
%                 *pScaleFactor;

        avArtifact(:,:)=avArtifact(:,:)...
            .*aScaleFactor.*pScaleFactor;
        zeroLine=max(avArtifact(avAfShift:end,electrodeRef+1));

        % plot artifact shape
        plot(avArtifact(avAfShift:end,electrodeRef+1)....
            -zeroLine,...
                        'color','red','XLimInclude','off');      
        
        %plot(originalData(...
        %        max(pulseTimes(1)-lOffset,1):min(pulseTimes(1)+pulseDuration,size(originalData))),'color','red');

    
    end
    
    % Some alignment
    if exist('yAxMin')&&exist('yAxMax')
        ylim([yAxMin yAxMax]);
    end
    
    % Add red line for pulse to plot
    title(strcat('Electrode',' ',num2str(electrodeRef)));
    xlabel('Samples');
    ylabel('pulse iterations (arbitrary units)');

    %ylim=get(gca,'ylim');
    line([lOffset;lOffset],ylim.',...
        'linewidth',1,...
        'color','red',...
        'YLimInclude','off');

    % Save plot
        position = get(gcf,'Position');
        set(gcf,'Color','w',...
        'PaperPositionMode', 'auto', ...
        'Units','in','Position',[position(1:2) 8 11 ],...
        'PaperPosition',[0.25 0.25 8 11])
        myStyle = hgexport('factorystyle');
        myStyle.Format = 'png';
        myStyle.Width = 8;
        myStyle.Height = 11;
        myStyle.Resolution = 600;
        myStyle.Units = 'inch';
        myStyle.FixedFontSize = 12;
    if savePlot
        hgexport(gcf,strcat(dataPath,'electrode_',num2str(electrodeRef)),myStyle,'Format','png')
    end
end

% xt=get(gca,'yticklabel');
% xt=str2num(xt);
% set(gca,'yticklabel',{'1','10','100','1000'});

if closeAllFigures;
    close all;
end

% Garbage collection
org.dt.matlab.utilities.JavaMemoryCleaner.clear(0);

