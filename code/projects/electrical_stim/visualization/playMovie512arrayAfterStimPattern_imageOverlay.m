function playMovie512arrayAfterStimPattern_imageOverlay(pathToAnalysisData, patternNo,imageData, electrodeCoords, varargin)
%% Show average recordings from all electrodes after a given stimulus from a single electrode as movies
% inputs:  pathToAnalysisData: a string that points to preprocessed data e.g.,'/Volumes/Analysis/2012-09-24-3/data008/';
%          patternNo: stimulation pattern
%          imageData: image (RGB or grayscale) to show underneat stimulation
%                   data 
%          electrodeCoords: electrode coordinates with respect to the image
%       
%    optional:   movieNo - only play a particular movie number or set of movies. omitting plays all movies available for a particular pattern.    
%                saveMovie - logical true or false
%                colorScale - controls clims on the plot, default [-20 10].
%                must be a vector of length 2 [lowerLim upperLim]
%                circleSize - default 350, controls the size of the marker
%                'electrodes'
%                saveInsideAxesOnly - allows saving only the inner axes,
%                excluding the surrounding whitespace. Must be used with
%                saveMovie set to true or it won't work. 
% Usage: playMovie512arrayAfterStimPattern_imageOverlay('/Volumes/Analysis/2012-09-24-3/data008/',9,'saveMovie',true)
% Lauren Grosberg 1/2016

% Get function arguments
p = inputParser;
p.addRequired('pathToAnalysisData', @ischar)
p.addRequired('patternNo', @isnumeric)

p.addParameter('movieNo', 0, @isnumeric) 
p.addParameter('movieIndex', 0, @isnumeric)
p.addParameter('saveMovie', false, @islogical) %default: don't save movie
p.addParameter('colorScale',[-20 10], @isnumeric); 
p.addParameter('circleSize', 350, @isnumeric);
p.addParameter('saveInsideAxesOnly',false,@islogical); 

p.parse(pathToAnalysisData, patternNo, varargin{:})
saveMovie = p.Results.saveMovie; 
movieNo = p.Results.movieNo;
movieIndex = p.Results.movieIndex;
colorScale = p.Results.colorScale; 
circleSize = p.Results.circleSize; 
saveInsideAxesOnly = p.Results.saveInsideAxesOnly; 

% Get electrode coordinates for the 512-electrode MEA from user input
positions = electrodeCoords';

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

if movieIndex
    mIndices = movieIndex;
end

dataTraces=NS_ReadPreprocessedData(pathToAnalysisData, '', 0, patternNo,...
    movieNos(1), 99999);

if saveMovie
    [movfile,movpath] = uiputfile('*.avi','Save Movie As');
    movieFileName = [movpath movfile];
    writerObj = VideoWriter(movieFileName);
    writerObj.FrameRate = 3;
    open(writerObj);
end
firstArtifact = mean(dataTraces,1);
f = figure; set(f,'Position',[100 465 845 445]);
set(f,'Color','white');
imshow(imageData); axis image; 
for movieIndex = mIndices
    
    dataTraces=NS_ReadPreprocessedData(pathToAnalysisData, '', 0, patternNo,...
        movieNos(movieIndex), 99999);
    % get stimulus amplitude
    [amps, stimChan, stimAmpVectors] = getStimAmps(pathToAnalysisData,...
        patternNo, movieNos(movieIndex));
    subtractionMatrix = repmat(firstArtifact,[size(dataTraces,1) 1]);
    
    for t = 1:40 %size(dataTraces,3)
        
        meanData = mean(dataTraces(:,:,t)-subtractionMatrix(:,:,t),1);
        colorsToPlot = meanData; 
        colorsToPlot((colorsToPlot<colorScale(1))) = colorScale(1); 
        colorsToPlot((colorsToPlot>colorScale(2))) = colorScale(2); 
        colorsToPlot = ceil(255*(colorsToPlot - colorScale(1))/diff(colorScale)); 
        cmap = cool; 
        colors = cmap(colorsToPlot+1,:); 
       
        if t<10
            meanData(stimChan) = 100*stimAmpVectors(:,2)';
        end
        if t==1
            hold on; 
            ah = scatter(positions(:,1),positions(:,2),10*abs(meanData),'filled');
            th = text(positions(508,1),positions(510,2),sprintf('%0.3f ms',t/20),'FontSize',16); 
        end
        hold on; delete(ah); 
        ah = scatter(positions(:,1),positions(:,2),8*abs(meanData)+0.1,colors); 
        ah.LineWidth = 4; 
        set(gca,'FontSize',16);
        
        title(sprintf('%s \npattern %0.0f; movie no. %0.0f; stimAmp %0.2f uA; t = %0.3f ms',pathToAnalysisData,patternNo,movieNos(movieIndex),amps(1),t/20));
        hold on; scatter(positions(stimChan,1),positions(stimChan,2),350,'black');
        delete(th); 
        th = text(positions(508,1),positions(510,2),sprintf('%0.3f ms',t/20),'FontSize',16); 
        
        if saveMovie
            if saveInsideAxesOnly
                title(''); %turn title off. 
                M = getframe(gca);
            else
                M = getframe(f);
            end
            writeVideo(writerObj,M);
        end
        pause(0.1); 
    end
    
end
if saveMovie
    close(writerObj);
end
end
