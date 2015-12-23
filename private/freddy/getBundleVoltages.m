function bundleMeans = getBundleVoltages(patternNos, display)
movieNo = 0;
% Load matrix containing the electrode numbers for the 512-electrode MEA
% temp = load('/Users/vision/Dropbox/Lab/Development/matlab-standard/private/freddy/512elecpositions.mat'); % Find a more general location for this or call a different text file.
temp = load('~/Development/matlab-standard/private/freddy/512elecpositions.mat'); % Find a more general location for this or call a different text file.

positions = temp.positions;

hexCoords = elecHexCoords();

bundleMeans = zeros(70,2, size(patternNos, 1));

% Find movie indices

for patternIndex = 1:size(patternNos, 2)
    movieNos = [];
    patternNo = patternNos(patternIndex);
    disp(['Testing pattern ' num2str(patternNo)]);
    pathToAnalysisData = '/Volumes/Analysis/2014-08-20-1/data003/';
    patternNoString = ['p' num2str(patternNos(patternIndex))];
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
    if display
        f = figure; set(f,'Position',[100 360 1000 550]);
        set(f,'Color','white');
    end
    
    for movieIndex = mIndices
        cla;
        dataTraces=NS_ReadPreprocessedData(pathToAnalysisData, '', 0, patternNo,...
            movieNos(movieIndex), 99999);
        % get stimulus amplitude
        [amps, stimChan, stimAmpVectors] = getStimAmps(pathToAnalysisData,...
            patternNo, movieNos(movieIndex));
        subtractionMatrix = repmat(firstArtifact,[size(dataTraces,1) 1]);
        
        
        bundle = zeros(16, 3);
        
        %for t = 1:30 %size(dataTraces,3)
        cla;
        [meanData, minIndices] = min(mean(dataTraces(:,:,10:30)-subtractionMatrix(:,:,10:30),1), [], 3);
        
        xpos = -945;
        ypos = -450;
        
        minID = electrodeIDForXY(positions, xpos, ypos);
        
        for xpos = -945:60:915
            testID = electrodeIDForXY(positions, xpos, ypos);
            if meanData(testID) < meanData(minID)
                minID = testID;
            end
        end
        
        bundle(1, 1) = meanData(minID);
        bundle(1, 2) = minID;
        bundle(1, 3) = minIndices(minID);
        
        xpos = positions(minID, 1);
        i = 2;
        for ypos = -390:60:450
            leftID = electrodeIDForXY(positions, xpos-30, ypos);
            %disp(['leftID: ' num2str(leftID) ' for x,y: ' num2str(xpos-30) ', ' num2str(ypos)]);
            rightID = electrodeIDForXY(positions, xpos+30, ypos);
            %disp(['rightID: ' num2str(rightID) ' for x,y: ' num2str(xpos+30) ', ' num2str(ypos)]);
            if leftID == 0 || (rightID && meanData(leftID) > meanData(rightID))
                xpos = xpos+30;
                bundle(i, 1) = meanData(rightID);
                bundle(i, 2) = rightID;
                bundle(i, 3) = minIndices(rightID);
            elseif (rightID == 0 && leftID) || (leftID && meanData(rightID) >= meanData(leftID))
                xpos = xpos-30;
                bundle(i, 1) = meanData(leftID);
                bundle(i, 2) = leftID;
                bundle(i, 3) = minIndices(leftID);
                
            end
            i = i + 1;
            
            
        end
        
        for j = 1:size(bundle, 1)
            if ismember(bundle(j, 2), stimChan) || ismember(bundle(j, 2), hexCoords)
                bundle(j, 1) = NaN;
            end
        end
        
        for j = 5:size(bundle, 1)
            %xpos = positions(bundle(j, 2), 1);
            %ypos = positions(bundle(j, 2), 2);
            
            atedge = true;
            backcount = 0;
            while (atedge && backcount < 4)
                xpos = positions(bundle(j-backcount, 2), 1);
                if mod(j-backcount,2) == 0
                    if xpos ~= 945
                        atedge = false;
                    end
                else
                    if xpos ~= 915
                        atedge = false;
                    end
                end
                backcount = backcount + 1;
            end
            if atedge
                disp(['At edge starting at row ' num2str(j)]);
                bundleTrimmed = zeros(j, 2);
                bundleTrimmed(1:j, :) = bundle(1:j, :);
                bundle = bundleTrimmed;
                break;
            end
     
        end
            
        
        bundleMean = mean(bundle(~isnan(bundle(:, 1))), 1);
        
        bundleMeans(movieIndex-1, 1, patternIndex) = bundleMean;
        bundleMeans(movieIndex-1, 2, patternIndex) = amps(1);
        
        if display
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