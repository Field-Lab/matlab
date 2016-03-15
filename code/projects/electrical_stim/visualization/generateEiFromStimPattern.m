function [rawData, amplitudes] = generateEiFromStimPattern(pathToAnalysisData, patternNo,varargin)
%% Show average recordings from all electrodes after a given stimulus as an EI
% inputs:  pathToAnalysisData: a string that points to preprocessed data e.g.,'/Volumes/Analysis/2012-09-24-3/data008/';
%          patternNo: numeric pattern number
%    optional:   movieNo - only play a particular movie number or set of movies. omitting plays all movies available for a particular pattern.    
%                saveImages - logical true or false
%                colorScale - controls clims on the plot, default [0 100].
%                must be a vector of length 2 [lowerLim upperLim]
%                circleSize - default 350, controls the size of the marker
%                'electrodes'
%                suppressPlots - option not to show plots, default 0
%                plotElecWaveforms
%                showElecNums logical true or false, default false, 
% Usage: generateEiFromStimPattern('/Volumes/Analysis/2012-09-24-3/data006/', 49,'movieNo',443)
% Lauren Grosberg 9/2014

% Get function arguments
p = inputParser;
p.addRequired('pathToAnalysisData', @ischar)
p.addRequired('patternNo', @isnumeric)

p.addParamValue('movieNo', 0, @isnumeric) 
p.addParamValue('saveImages', false, @islogical) %default: don't save movie
p.addParamValue('colorScale',[0 1000], @isnumeric); 
p.addParamValue('circleSize', 350, @isnumeric);
p.addParamValue('showElecNums',false, @islogical);
p.addParamValue('plotElecWaveforms',[]);
p.addParamValue('suppressPlots',false,@islogical);

p.parse(pathToAnalysisData, patternNo, varargin{:})
saveImages = p.Results.saveImages; 
movieNo = p.Results.movieNo;
colorScale = p.Results.colorScale; 
circleSize = p.Results.circleSize; 
showElecNums = p.Results.showElecNums; 
plotElecWaveforms = p.Results.plotElecWaveforms; 
suppressPlots = p.Results.suppressPlots; 

% Load matrix containing the electrode numbers for the 512-electrode MEA
currentFunctionPath = mfilename('fullpath');
temp = load(fullfile(currentFunctionPath, '../../resources/arrayPositions512.mat'));
positions = temp.positions;

if ~strcmp(pathToAnalysisData(end),filesep)
    pathToAnalysisData = [pathToAnalysisData filesep];
end

% Find movie indices
movieNos = findMovieNos(pathToAnalysisData,patternNo);

if movieNo
    mIndices = find(movieNos == movieNo);
else
    mIndices = 2:size(movieNos,2);
end

dataTraces=NS_ReadPreprocessedData(pathToAnalysisData, '', 0, patternNo,...
    movieNos(1), 99999);
firstArtifact = mean(dataTraces,1);
if ~suppressPlots
    f = figure; set(f,'Position',[100 465 845 445]);
    set(f,'Color','white');
end
rawData = zeros(size(dataTraces,1),size(dataTraces,2),size(dataTraces,3),length(mIndices));
idx = 0; % index for saving rawData output
for movieIndex = mIndices
    idx = idx + 1; 
    if ~suppressPlots; cla; end %or clear specific axis. 
    dataTraces=NS_ReadPreprocessedData(pathToAnalysisData, '', 0, patternNo,...
        movieNos(movieIndex), 99999);
    try
        rawData(:,:,:,idx) = dataTraces(1:size(rawData,1),1:size(rawData,2),1:size(rawData,3));
    catch
        rawData(end,:,:,:) = []; %Sometimes certain movies have fewer repetitions than the first one
        rawData(:,:,:,idx) = dataTraces(1:size(rawData,1),1:size(rawData,2),1:size(rawData,3));
    end
    % get stimulus amplitude
    [amps, stimChan, ~] = getStimAmps(pathToAnalysisData,...
        patternNo, movieNos(movieIndex));
    subtractionMatrix = repmat(firstArtifact,[size(dataTraces,1) 1]);
    modData = dataTraces - subtractionMatrix;
%     modData = modData(:,:,8:end); %cut off the early time points
%     spiketrials = [2 3 4 6 10 11 14 15 16 17 18 19 21 23]; %banana delete this line!!!!
%     modData = modData(spiketrials,:,:); 
    trialAvg = squeeze(mean(modData)); 
%     keyboard; 
    amplitudes = max(trialAvg,[],2) - min(trialAvg,[],2);
    scatterData = amplitudes;
    if ~suppressPlots
        %{
        scatter(positions(:,1),positions(:,2),circleSize,scatterData,'filled');
        axis off; axis image; c=colorbar;
        caxis(colorScale);
        ylabel(c,'  \muV','rot',0); set(gca,'FontSize',16);
        %}
        scatterData = abs(amplitudes) + 0.1;
        
        f0 = figure; set(f0,'Position',[50 465 845 445]);
        set(f0,'Color','white');
        [~,eiM] = ei2matrix(scatterData);
        imagesc(eiM,[0 max(eiM(:))/2]);
        axis image; axis off; colorbar; colormap hot;
        title(sprintf('%s \npattern %0.0f; movie no. %0.0f; stimAmp %0.2f uA',pathToAnalysisData,patternNo,movieNos(movieIndex),amps(1)));
        
       % {
        f2 = figure; set(f2,'Position',[100 465 845 445]);
        set(f2,'Color','white');
        scatter(positions(:,1),positions(:,2),scatterData,'filled','green','ButtonDownFcn',{@plotWaveform,positions});
        axis image; axis off;
        title(sprintf('%s \npattern %0.0f; movie no. %0.0f; stimAmp %0.2f uA',pathToAnalysisData,patternNo,movieNos(movieIndex),amps(1)));
        hold on; scatter(positions(stimChan,1),positions(stimChan,2),350,'black');
        text(positions(stimChan,1),positions(stimChan,2),'stimulating electrode');
       % }
        
        if showElecNums % Display electrode numbers
            for e = 1:512
                text(positions(e,1),positions(e,2),num2str(e),'HorizontalAlignment','center');
            end
        end
        if plotElecWaveforms
            f3 = figure; set(f3, 'Color', 'white');
            numPlots = length(plotElecWaveforms);
            for nP = 1:numPlots
                figure(f3);
                if numPlots <= 5
                    subplot(1,numPlots,nP);
                elseif numPlots <= 10
                    subplot(2,ceil(numPlots/2),nP);
                elseif numPlots <= 15
                    subplot(3,ceil(numPlots/3),nP);
                end
                plot(linspace(0,size(trialAvg,2)/20,size(trialAvg,2)), trialAvg(plotElecWaveforms(nP),:));
                title(sprintf('electrode %0.0f\n %0.0f', plotElecWaveforms(nP), nP));
                xlabel('ms');  xlim([0 1.5]); ylim([-100 100]);
                
                figure(f2);
                text(positions(plotElecWaveforms(nP),1), ...
                    positions(plotElecWaveforms(nP),2), num2str(nP), ...
                    'HorizontalAlignment','center');
            end
        end %plot waveforms
    end %suppressPlots
end

if saveImages
    [imfile,impath] = uiputfile('*.*','Save Images As');
    imageFileName = [impath imfile];
    saveas(f,imageFileName,'epsc'); saveas(f,imageFileName,'jpeg');
    saveas(f2,imageFileName,'epsc'); saveas(f2,imageFileName,'jpeg');
end

    function plotWaveform(src,~,positions)
        axesHandle = get(src,'Parent');
        a = get(axesHandle,'CurrentPoint');
        x = a(1,1);
        y = a(1,2);
        diffX = positions(:,1) - x;
        diffY = positions(:,2) - y;
        d = hypot(diffX,diffY);
        electrode = find(d == min(d));
        figure;
        plot(linspace(0,size(trialAvg,2)/20,size(trialAvg,2)), trialAvg(electrode,:));
        title(sprintf('electrode %0.0f', electrode));
        xlabel('ms');  xlim([0 1.5]); ylim([-100 100]);
        figure(f2);
        text(positions(electrode,1), ...
            positions(electrode,2), num2str(electrode), ...
            'HorizontalAlignment','center');
%         keyboard;
%         set(src,'CData',[0 0 1]); 
    end
end