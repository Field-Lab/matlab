function bundleMeans = getBundleVoltagesAStar(pathToAnalysisData, patternNos, display, varargin)
% Calculate the average bundle voltage for each movie of each pattern in 
% patternNos, using the A* pathfinding algorithm to locate the bundle
% Optional inputs

nbin = length(varargin);
if mod(nbin,2)==1
    err = MException('MATLAB:InvArgIn','Unexpected number of arguments');
    throw(err);
end

movieNo = 0;
% Load matrix containing the electrode numbers for the 512-electrode MEA
temp = load([matlab_code_path() ...
    'code/projects/electrical_stim/resources/arrayPositions512.mat']); 
positions = temp.positions;

% Set up default parameters. 
exclude = 0;
movieIndex = []; 
% Read the optional input arguments
for j=1:(nbin/2)
    if ~ischar(varargin{j*2-1})
        err = MException('MATLAB:InvArgIn',...
            'Unexpected additional property');
        throw(err);
    end
    
    switch lower(varargin{j*2-1})
        case 'exclude'
            exclude = varargin{j*2};
        case 'movieindex'
            movieIndex = varargin{j*2};

        otherwise
            err = MException('MATLAB:InvArgIn',...
                'Unknown parameter specified');
            throw(err);
    end
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
    
    patternNo = patternNos(patternIndex);
    disp(['Testing pattern ' num2str(patternNo)]);
 
    movieNos = findMovieNos(pathToAnalysisData,patternNo);
    if isempty(movieIndex)
        mIndices = 2:size(movieNos,2);
    else
        mIndices = movieIndex;
    end
    dataTraces=NS_ReadPreprocessedData(pathToAnalysisData, '', 0, patternNo,...
        movieNos(1), 99999);
    
    firstArtifact = mean(dataTraces,1);
    if display
        f = figure; set(f,'Position',[100 360 1000 550]);
        set(f,'Color','white');
    end
    
    ignore = [];
    
    for movieIndex = mIndices
        dataTraces=NS_ReadPreprocessedData(pathToAnalysisData, '', 0, patternNo,...
            movieNos(movieIndex), 99999);
        % get stimulus amplitude
        [amps, stimChan, stimAmpVectors] = getStimAmps(pathToAnalysisData,...
            patternNo, movieNos(movieIndex));
        subtractionMatrix = repmat(firstArtifact,[size(dataTraces,1) 1]);
        
        if isempty(ignore)
            ignore = hexNeighborsFast(stimChan, hexArray, hexCoords);
        end
        
        %for t = 1:30 %size(dataTraces,3)
        [meanData, minIndices] = min(mean(dataTraces(:,:,10:30)-subtractionMatrix(:,:,10:30),1), [], 3);
        
        axonPath = findAxonBundlePathFast(meanData, borderElecs, hexCoords, 0, 0, exclude);
        bundle = zeros(size(axonPath, 2), 3);
        bundle(:, 1) = meanData(axonPath);
        bundle(:, 2) = axonPath;
        bundle(:, 3) = minIndices(axonPath);
                
        for j = 1:size(bundle, 1)
            if ismember(bundle(j, 2), ignore)
                bundle(j, 1) = NaN;
            end
        end
        
         
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
        if 0 % bypass for now.
        if amps(1) < -3.5
            figure; scatter(positions(:,1),positions(:,2),abs(meanData),'filled')
            hold on; scatter(positions(axonPath,1),positions(axonPath,2),100,'k')
            hold on; scatter(positions(ignore,1),positions(ignore,2),100,'r');
            axis image; axis off; title('black circles were averaged, red ignored');
            title(sprintf(['%s \npattern %0.0f; movie no. %0.0f; stimAmp %0.2f uA',...
                '\nblack circles were averaged, red ignored'],pathToAnalysisData,...
                patternNo,movieNos(movieIndex),amps(1)), 'Color', 'black');  
        end
        end
        
    end
end
end