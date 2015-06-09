function [threshold,completeFit, erfErr] = fitToErf(elecResp,plotResponseCurves) 

if isfield(elecResp, 'stimInfo') %if elecResp is loaded up
    nMovies = length(elecResp.stimInfo.movieNos);
      
    data = zeros(2, nMovies);
    data(1,:) = elecResp.stimInfo.stimAmps;
    data(2,:) = elecResp.analysis.successRates;
    data(3,:) = elecResp.stimInfo.nPulses;
    lockedAmps = elecResp.analysis.finalized;
    for i = length(elecResp.stimInfo.stimAmps): -1: 1
        if isempty(elecResp.analysis.type{i})
            data(:,i) = [];
            lockedAmps(i) = [];
        end
    end
else
    warnH = warndlg('need to load up an elecResp file first!');
end

% linear-based
data(1,:) = abs(data(1,:));
[erfParams, completeFit, erfErr] = erfFitter(data, 2, -1, 'makePlot', plotResponseCurves, 'lockedAmps', lockedAmps);
threshold = -erfParams(2)/erfParams(1);

%standard deviations: don't redo if analysis hasn't changed

% threshStd = bootstrapThresh(elecResp, main.bootstrapReps);
% bootstrapReps = main.bootstrapReps;

end

