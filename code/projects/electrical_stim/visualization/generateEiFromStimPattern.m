function generateEiFromStimPattern(pathToAnalysisData, patternNo,varargin)
%% Show average recordings from all electrodes after a given stimulus as an EI
% inputs:  pathToAnalysisData: a string that points to preprocessed data e.g.,'/Volumes/Analysis/2012-09-24-3/data008/';
%          patternNo: numeric pattern number
%    optional:   movieNo - only play a particular movie number or set of movies. omitting plays all movies available for a particular pattern.    
%                saveImages - logical true or false
%                colorScale - controls clims on the plot, default [0 100].
%                must be a vector of length 2 [lowerLim upperLim]
%                circleSize - default 350, controls the size of the marker
%                'electrodes'
% Usage: generateEiFromStimPattern('/Volumes/Analysis/2012-09-24-3/data006/', 49,'movieNo',443)
% Lauren Grosberg 9/2014

% Get function arguments
p = inputParser;
p.addRequired('pathToAnalysisData', @ischar)
p.addRequired('patternNo', @isnumeric)

p.addParamValue('movieNo', 0, @isnumeric) 
p.addParamValue('saveImages', false, @islogical) %default: don't save movie
p.addParamValue('colorScale',[0 100], @isnumeric); 
p.addParamValue('circleSize', 350, @isnumeric);
p.addParamValue('showElecNums',false, @islogical);
p.addParamValue('plotElecWaveforms',[]);

p.parse(pathToAnalysisData, patternNo, varargin{:})
saveImages = p.Results.saveImages; 
movieNo = p.Results.movieNo;
colorScale = p.Results.colorScale; 
circleSize = p.Results.circleSize; 
showElecNums = p.Results.showElecNums; 
plotElecWaveforms = p.Results.plotElecWaveforms; 

% Load matrix containing the electrode numbers for the 512-electrode MEA
temp = load('/Users/grosberg/matlab/arrayPositions512.mat'); % Find a more general location for this or call a different text file.
positions = temp.positions;

if ~strcmp(pathToAnalysisData(end),filesep)
    pathToAnalysisData = [pathToAnalysisData filesep];
end

% Find movie indices
movieNos = [];
patternNoString = ['p' num2str(patternNo)];
files = dir([pathToAnalysisData patternNoString]);

for i = 1:length(files)
    if strfind(files(i).name, patternNoString) == 1
        mIndices = strfind(files(i).name, 'm');
        movieNos = [movieNos str2double(files(i).name(mIndices(end)+1:end))]; %#ok<AGROW>
    end
end

movieNos = sort(movieNos);
if movieNo
    mIndices = find(movieNos == movieNo);
else
    mIndices = 2:size(movieNos,2);
end

dataTraces=NS_ReadPreprocessedData(pathToAnalysisData, '', 0, patternNo,...
    movieNos(1), 99999);
firstArtifact = mean(dataTraces,1);
f = figure; set(f,'Position',[100 465 845 445]);
set(f,'Color','white');
for movieIndex = mIndices
    cla;
    dataTraces=NS_ReadPreprocessedData(pathToAnalysisData, '', 0, patternNo,...
        movieNos(movieIndex), 99999);
    % get stimulus amplitude
    [amps, stimChan, ~] = getStimAmps(pathToAnalysisData,...
        patternNo, movieNos(movieIndex));
    subtractionMatrix = repmat(firstArtifact,[size(dataTraces,1) 1]);
    modData = dataTraces - subtractionMatrix;
%     modData = modData(:,:,8:end); %cut off the early time points
    trialAvg = squeeze(mean(modData)); 
    amplitudes = max(trialAvg,[],2) - min(trialAvg,[],2);
    scatterData = amplitudes;
    scatter(positions(:,1),positions(:,2),circleSize,scatterData,'filled');
    axis off; axis image; c=colorbar;
    caxis(colorScale);
    ylabel(c,'  \muV','rot',0); set(gca,'FontSize',16);
    
    scatterData = abs(amplitudes) + 0.1;
    
    f0 = figure; set(f0,'Position',[50 465 845 445]);
    set(f0,'Color','white');
    [~,eiM] = ei2matrix(scatterData); 
    imagesc(eiM,[0 max(eiM(:))/2]); 
    axis image; axis off; colorbar; colormap hot;
%     title(['EI for n ' num2str(neuron)]); axis off;
    title(sprintf('%s \npattern %0.0f; movie no. %0.0f; stimAmp %0.2f uA',pathToAnalysisData,patternNo,movieNos(movieIndex),amps(1)));
  
    
     
    f2 = figure; set(f2,'Position',[100 465 845 445]);
    set(f2,'Color','white');
    scatter(positions(:,1),positions(:,2),scatterData,'filled','green');
    axis image; axis off;
    title(sprintf('%s \npattern %0.0f; movie no. %0.0f; stimAmp %0.2f uA',pathToAnalysisData,patternNo,movieNos(movieIndex),amps(1)));
    hold on; scatter(positions(stimChan,1),positions(stimChan,2),350,'black');
    text(positions(stimChan,1),positions(stimChan,2),'stimulating electrode');
    
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
    end
end

if saveImages
    [imfile,impath] = uiputfile('*.*','Save Images As');
    imageFileName = [impath imfile];
    saveas(f,imageFileName,'epsc'); saveas(f,imageFileName,'jpeg');
    saveas(f2,imageFileName,'epsc'); saveas(f2,imageFileName,'jpeg');
end
end