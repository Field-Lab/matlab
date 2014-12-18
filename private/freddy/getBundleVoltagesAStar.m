function bundleMeans = getBundleVoltagesAStar(path, patternNos, display, exclude)
%Calculate the average bundle voltage for each movie of each pattern in patternNos, using the A* pathfinding algorithm to locate the bundle.
movieNo = 0;
% Load matrix containing the electrode numbers for the 512-electrode MEA
temp = load('~/git_code/matlab/private/freddy/512elecpositions.mat'); % Find a more general location for this or call a different text file.
positions = temp.positions;

if nargin < 4
    exclude = 0;
end

topBorder = [249:256 261:8:381 385:392];
rightBorder = 392:8:512;
bottomBorder = [505:512 8:8:136 129:135];
leftBorder = 129:8:249;
borderElecs = unique([topBorder rightBorder bottomBorder leftBorder]);

hexCoords = elecHexCoords();
hexArray = elecHexArray(hexCoords);

bundleMeans = zeros(70,3, size(patternNos, 1));

% Find movie indices

for patternIndex = 1:size(patternNos, 2)
    movieNos = [];
    patternNo = patternNos(patternIndex);
    disp(['Testing pattern ' num2str(patternNo)]);
    pathToAnalysisData = path; %'/Volumes/Analysis/2012-09-24-3/data008/';
    patternNoString = ['p' num2str(patternNos(patternIndex))];
    files = dir([pathToAnalysisData patternNoString]);
 
    for i = 1:length(files)
        if strfind(files(i).name, patternNoString) == 1
            mIndices = strfind(files(i).name, 'm');
            movieNos = [movieNos str2double(files(i).name(mIndices(end)+1:end))]; %#ok<AGROW>
        end
    end
    movieNos = sort(movieNos);
    mIndices = 2:size(movieNos,2);
    dataTraces=NS_ReadPreprocessedData(pathToAnalysisData, '', 0, patternNo,...
        movieNos(1), 99999);
    
    firstArtifact = mean(dataTraces,1);
    if display
        f = figure; set(f,'Position',[100 360 1000 550]);
        set(f,'Color','white');
    end
    
    for movieIndex = mIndices
        dataTraces=NS_ReadPreprocessedData(pathToAnalysisData, '', 0, patternNo,...
            movieNos(movieIndex), 99999);
        % get stimulus amplitude
        [amps, stimChan, stimAmpVectors] = getStimAmps(pathToAnalysisData,...
            patternNo, movieNos(movieIndex));
        subtractionMatrix = repmat(firstArtifact,[size(dataTraces,1) 1]);
                
        %for t = 1:30 %size(dataTraces,3)
        [meanData, minIndices] = min(mean(dataTraces(:,:,10:30)-subtractionMatrix(:,:,10:30),1), [], 3);
        
        path = findAxonBundlePathFast(meanData, borderElecs, hexCoords, 0, 0, exclude);
        bundle = zeros(size(path, 2), 3);
        bundle(:, 1) = meanData(path);
        bundle(:, 2) = path;
        bundle(:, 3) = minIndices(path);
        
        for j = 1:size(bundle, 1)
            if ismember(bundle(j, 2), stimChan) || ismember(bundle(j, 2), hexNeighborsFast(stimChan, hexArray, hexCoords))
                bundle(j, 1) = NaN;
            end
        end
        
%         for j = 5:size(bundle, 1)            
%             atedge = true;
%             backcount = 0;
%             while (atedge && backcount < 4)
%                 xpos = positions(bundle(j-backcount, 2), 1);
%                 if mod(j-backcount,2) == 0
%                     if xpos ~= 945
%                         atedge = false;
%                     end
%                 else
%                     if xpos ~= 915
%                         atedge = false;
%                     end
%                 end
%                 backcount = backcount + 1;
%             end
%             if atedge
%                 disp(['At edge starting at row ' num2str(j)]);
%                 bundleTrimmed = zeros(j, 2);
%                 bundleTrimmed(1:j, :) = bundle(1:j, :);
%                 bundle = bundleTrimmed;
%                 break;
%             end
%      
%         end
            
        
        bundleMean = mean(bundle(~isnan(bundle(:, 1))), 1);
        
        bundleMeans(movieIndex-1, 1, patternIndex) = bundleMean;
        bundleMeans(movieIndex-1, 2, patternIndex) = amps(1);
        bundleMeans(movieIndex-1, 3, patternIndex) = movieNos(movieIndex);
        
        if display
            cla;
            scatter(positions(:,1),positions(:,2),350,meanData,'filled');
            axis off; axis image; colorbar;
            caxis([-50 -10]);
            
            title(sprintf('%s \npattern %0.0f; movie no. %0.0f; stimAmp %0.2f uA',pathToAnalysisData,patternNo,movieNos(movieIndex),amps(1)), 'Color', 'black');
            hold on; scatter(positions(bundle(:, 2),1),positions(bundle(:, 2),2),350,[0.5 0.5 0.5], 'filled');
            %         text(positions(stimChan,1),positions(stimChan,2),'stimulating electrode')
            pause(0.001);
        end
        %end
        
    end
end
end