function playMovie512arrayAfterStimPattern(pathToAnalysisData, patternNo,saveMovie)
%% Show average recordings from all electrodes after a given stimulus from a single electrode as movies
% inputs:  pathToAnalysisData: a string that points to preprocessed data e.g.,'/Volumes/Analysis/2012-09-24-3/data008/';
%          patternNo: 
%          saveMovie


% Load matrix containing the electrode numbers for the 512-electrode MEA
fname = fullfile(fileparts(mfilename('fullpath')),'../../resources/array_matrix_id510'); 
h = load(fname); 
electrodeMatrix = h.array_matrix_id510; clear h; 

if ~strcmp(pathToAnalysisData(end),filesep)
    pathToAnalysisData = [pathToAnalysisData filesep];
end

data_path = pathToAnalysisData;
% Find movie indices
movieNos = [];
patternNoString = ['p' num2str(patternNo)];
files = dir([data_path patternNoString]);

for i = 1:length(files)
    if strfind(files(i).name, patternNoString) == 1
        mIndices = strfind(files(i).name, 'm');
        movieNos = [movieNos str2double(files(i).name(mIndices(end)+1:end))]; %#ok<AGROW>
    end
end
movieNos = sort(movieNos);

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
f = figure; set(f,'Position',[100 360 1000 550]);
for movieIndex = 1:size(movieNos,2)
    cla;
    dataTraces=NS_ReadPreprocessedData(pathToAnalysisData, '', 0, patternNo,...
        movieNos(movieIndex), 99999);
    % get stimulus amplitude
    [amps, electrodes, stimAmpVectors] = getStimAmps(pathToAnalysisData,...
        patternNo, movieNos(movieIndex));
    subtractionMatrix = repmat(firstArtifact,[size(dataTraces,1) 1]);
    patternPath = [pathToAnalysisData 'pattern_files/'];
    temp = dir([patternPath 'pattern' num2str(patternNo) '_m*.mat']);
    temp = load([patternPath temp(1).name]);
    if isfield(temp,'pattern')
        for ii = 1:size(temp.pattern,2)
            stimChan(1,ii) = temp.pattern(ii).channel;
        end
    elseif isfield(temp,'Pattern')
        for ii = 1:size(temp.Pattern,2)
            stimChan(1,ii) = temp.Pattern(ii).channel;
        end
    end
    y = zeros(1,size(stimChan,2));  x = zeros(1,size(stimChan,2));
    for ii = 1:size(stimChan,2)
        [y(ii),x(ii)] = find(electrodeMatrix == stimChan(ii));
    end
    for t = 1:40 %size(dataTraces,3)
        cla;
        %             subplot(3,4,t-7);
        meanData = mean(dataTraces(:,:,t)-subtractionMatrix(:,:,t),1);
        if t<5
            meanData(stimChan) = 100;
        end
        [~,EIm_view]   = ei2matrix(meanData);
        imagesc(EIm_view,[-20 10]); axis image;  axis off; colorbar;
        
        title(sprintf('%s \npattern %0.0f; movie no. %0.0f; stimAmp %0.2f uA; t = %0.3f ms',pathToAnalysisData,patternNo,movieNos(movieIndex),amps(1),t/20));
        if round(x/2) == x/2
            text(x*2+2,y*2,'stimulating electrode')
        else
            text(x*2+1,y*2,'stimulating electrode')
        end
        if saveMovie
            M = getframe(f);
            writeVideo(writerObj,M);
        end
        pause(0.05);
    end
    
end
if saveMovie
    close(writerObj);
end
end


