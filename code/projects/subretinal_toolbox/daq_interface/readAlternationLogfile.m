function [ gratingsMap ] = readAlternationLogfile(alternationLogfilePath)
% READALTERNATIONLOGFILE 
%
% Outputs a map such that the set of keys is the set of grating sizes, 
% and the set of values is 2xn where the first line is the set of phases
% done for the particular grating, and the second line is the experiment ID
% corresponding to that particular size/phase combination.

gratingsMap = containers.Map('KeyType', 'int32', 'ValueType', 'any');

gratingSizeCol = 3;
gratingPhaseCol = 7;
logfileData = dlmread(alternationLogfilePath, ',', 2, 0);
gratData = int32(logfileData(:,gratingSizeCol));
phaseData = logfileData(:,gratingPhaseCol);
nLines = size(logfileData, 1);

cExpID = 0;
cGratSize = -1;
cGratPhase = -1;
for kk=1:nLines
    newGratSize = gratData(kk);
    newGratPhase = phaseData(kk);
    
    if (cGratPhase ~= newGratPhase) || (cGratSize ~= newGratSize)
        cExpID = cExpID + 1;
        if isKey(gratingsMap, newGratSize)
            gratingsMap(newGratSize) = [gratingsMap(newGratSize) [newGratPhase; cExpID]];
        else
            gratingsMap(newGratSize) = [newGratPhase; cExpID];
        end
    end
    
    cGratSize = newGratSize;
    cGratPhase = newGratPhase;
end

end

