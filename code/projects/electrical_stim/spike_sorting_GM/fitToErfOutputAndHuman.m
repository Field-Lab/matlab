function [thresholdHum thresholdAlg] = fitToErfOutputAndHuman(Output) 
%The same as fitToErf but now we make a comparison with the spike sorting
%algorithm using information provided by the Output structure (no need to
%specify elecResp)
%Gonzalo Mena 06/15
pathToAnalysisData = Output.path;
patternNo = Output.stimInfo.patternNo;
neuronIds = Output.neuronInfo.neuronIds;
J = Output.tracesInfo.J;
I = Output.tracesInfo.I;
for n=1:length(neuronIds)
    tempPath = [pathToAnalysisData 'elecResp_n' num2str(neuronIds(n)) '_p' num2str(patternNo) '.mat']; 
    load(tempPath)
    

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
[erfParams, completeFit, erfErr] = erfFitter(data, 2, -1, 'makePlot', 0, 'lockedAmps', lockedAmps);
thresholdHum(n) = -erfParams(2)/erfParams(1);

data2(1,:)=abs(Output.stimInfo.listAmps(:,1))';

for j=1:J
    data2(2,j)=nansum(Output.spikes{n}(j,:))/I(j);
end

data2(3,:)=I;



[erfParams, completeFit, erfErr] = erfFitter(data2, 2, -1, 'makePlot', 0, 'lockedAmps', lockedAmps);
thresholdAlg = -erfParams(2)/erfParams(1);

%standard deviations: don't redo if analysis hasn't changed

% threshStd = bootstrapThresh(elecResp, main.bootstrapReps);
% bootstrapReps = main.bootstrapReps;

end

