function [threshold,completeFit, erfErr]  = getThresholdFromAuto(elecRespAuto)
% GETTHRESHOLDFROMAUTO(..) returns the 50% activation threshold from an
% elecRespAuto file.
% 02/2016

nMovies = size(elecRespAuto.LogisticReg,2);
data = zeros(2, nMovies);
data(1,:) = elecRespAuto.stimInfo.listAmps;
data(2,:) = elecRespAuto.LogisticReg;
data(3,:) = elecRespAuto.tracesInfo.I;

% linear-based
data(1,:) = abs(data(1,:));
[erfParams, completeFit, erfErr] = erfFitter(data, 2, -1, 'makePlot', 0);
threshold = -erfParams(2)/erfParams(1);
