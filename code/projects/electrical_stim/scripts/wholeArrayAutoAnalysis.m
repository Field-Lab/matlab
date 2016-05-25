function wholeArrayAutoAnalysis(eiFilePath, pathToAnalysisData,varargin)
% Uses Gonzalo's spike sorting code to automatically analyze all the
% neurons in a preparation. 
% Inputs:       eiFilePath (e.g. /Volumes/Analysis/2015-10-06-3/data000/data000.ei)
%               pathToAnalysisData (e.g. /Volumes/Analysis/2015-10-06-3/data001-data002/)
%   optional    threshold
%               sortingThreshold
%               startPattern
%               endPattern
% LG 12/2015

% Set default parameters
startPattern = 1; 
endPattern = 512; 
threshold = 10; 
sortingThreshold = 45; % Difficult to analyze cells underneath this threshold

nbin = length(varargin);
if mod(nbin,2)==1
    err = MException('MATLAB:InvArgIn','Unexpected number of arguments');
    throw(err);
end

% Read the optional input arguments
for j=1:(nbin/2)
    if ~ischar(varargin{j*2-1})
        err = MException('MATLAB:InvArgIn',...
            'Unexpected additional property');
        throw(err);
    end
    
    switch lower(varargin{j*2-1})
        case 'startpattern'
            startPattern = varargin{j*2};
        case 'threshold'
            threshold = varargin{j*2};
        case 'sortingthreshold'
            sortingThreshold = varargin{j*2};
        case 'endpattern'
            endPattern = varargin{j*2};
        otherwise
            err = MException('MATLAB:InvArgIn',...
                'Unknown parameter specified');
            throw(err);
    end
end

% Find pattern-cell pairs to analyze.
disp('Finding pattern-cell pairs to analyze...'); 
[eiM,nIdList] = convertEiFileToMatrix(eiFilePath);
cellsToAnalyze = cell(512,1); 
for e = 1:512
    electrodeNo = e;
    cellIDs = [];
    for n = 1:1:size(nIdList,1)
        tempID  = nIdList(n);          % gets the id number of the ith neuron
        ei = squeeze(eiM(:,:,n))';
%         eiAmps = max(ei)-min(ei);
        eiMin = abs(min(ei)); 
        if eiMin(electrodeNo) > threshold && (max(eiMin(:)) > sortingThreshold)
            cellIDs = cat(1,cellIDs,tempID);
        end
%         if eiAmps(electrodeNo) > threshold
%             cellIDs = cat(1,cellIDs,tempID);
%         end
    end
    cellsToAnalyze{e} = cellIDs; 
    fprintf('*');
    if mod(e,80) == 0
        fprintf('electrode %0.0f \n',e); 
    end
end
fprintf('\n'); 


% %% %% % % %% %% % % %% %% % % %% %% % % %% %% % % %% %% % % %% %% % % %% 
% %% %% % Calculate the total number of pattern-cell pairs to analyze % %% 
% %% %% % % %% %% % % %% %% % % %% %% % % %% %% % % %% %% % % %% %% % % %% 
totalNo = 0; 
for e = startPattern:endPattern
    totalNo = length(cellsToAnalyze{e}) + totalNo; 
end
fprintf('Total number of pattern-cell pairs is %0.0f\n',totalNo); 

% %% %% % % %% %% % % %% %% % % %% %% % % %% %% % % %% %% % % %% %% % % %% 
% %% %% % Update to reflect the latest pattern-cell pair
% %% %% % % %% %% % % %% %% % % %% %% % % %% %% % % %% %% % % %% %% % % %% 
currentPair = 0; 
if startPattern > 1
    for e = 1:startPattern-1
        currentPair = length(cellsToAnalyze{e}) + currentPair;
    end
end

% %% %% % % %% %% % % %% %% % % %% %% % % %% %% % % %% %% % % %% %% % % %% 
% %% %% % Begin auto analysis % %% %% % 
% %% %% % % %% %% % % %% %% % % %% %% % % %% %% % % %% %% % % %% %% % % %% 
disp('Beginning automatic analysis. Come back in a while'); 
try
    for e = startPattern:endPattern
        cellIDs = cellsToAnalyze{e} ;
        patternNo = e;
        disp( '**************************************************');
        disp(['****** analyzing pattern no. ' num2str(e) ' *****************']);
        disp( '**************************************************');
        
        for n = 1:length(cellIDs)
            currentPair = 1+currentPair;
            nId = cellIDs(n);
            elecStimSorting(pathToAnalysisData,eiFilePath,patternNo,nId);
            disp(['finished analyzing neuron ' num2str(nId) ' for pattern '...
                num2str(patternNo) ' : ' num2str(currentPair) ' of ' num2str(totalNo)]);
        end
    end
catch ME
    disp(['Final electrode analyzed was ' num2str(e)'])
    rethrow(ME)
end