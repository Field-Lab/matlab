function elecResp = checkForUnfinishedAnalysis(elecResp, bootstrapReps, varargin)

% checks to make sure that all higher-order analysis values are based are current and resolves some
% differences in elecResp structure between the current and older versions of the analysis code
%
% checks for incomplete analysis:
%   (if checkForAnalysisGaps == true)
%   warns if there are any stimulus amplitudes that haven't been analyzed falling in
%   between other amplitudes that have been analyzed (i.e. not before beginning or after end of
%   analyzed curve)
%
%   warns if there are any analyzed, but 'unlocked', stimulus amplitudes
%
% checks for current higher-order analysis values:
%   makes sure that values for the following exist, and recalculates them if
%   elecResp.analysis.erfCurrent == 0 (signifying spike-sorting has been done since last
%   calculation of values)
%
%   values checked for: (elecResp.analysis.X)
%     erfParams (doesn't check for existence; recalculates only if erfCurrent == 0 --> also updates
%       threshold, erfErr)
%     threshStd (checks for existence-->also updates details.bootstrapReps)
%
%     --beging phased out--
%     logErfParams (checks for existence-->also updates, logThresh, logErfErr)
%     logThreshStd (checks for existence-->also updates details.logBootstrapReps)
%
%
%   ***note: prior to 2009-11-23, elecResp.analysis.erfCurrent didn't reflect the accuracy of all of
%   these higher-order values, so this function should be run with recalcAll == 1 before trusting
%   .erfCurrent for elecResp structures created before this date
%
%   ***note: erfCurrent flag doesn't reflect accuracy of constrainedSlope-related values
%
% updates to elecResp structure due to updates in code:
%
%   elecResp.stimInfo.pElec should specify primary electrode for multi-electrode stimulation patterns
%   elecResp.stimInfo.stimAmps should reflect amplitude on primary electrode if specified
%   elecResp.stimInfo.nPulses should be set separately for each movie (should be a vector of values)
%   elecResp.stimInfo.successRates should be correct (was an error in old code version)
%
% updated Dec. 2010: changed limit on number of traces loaded by NS_readPreprocessedData from 100 to
% 99999
%
% updated June 2011: changed default of 'keepLogBased' from 1 to 0
% updated July 2011: added parameter to suppress fitting with
% gaussCdfFitter in addition to erfFitter, and gaussCdfFitter parameters
% removed from elecResp if present



p = inputParser;

p.addRequired('elecResp', @isstruct)
p.addRequired('bootstrapReps', @isnumeric)

p.addParamValue('recalcAll', 0, @(x)any(x==[0 1]))
%p.addParamValue('keepLogBased', 0, @(x)any(x==[0 1]))
p.addParamValue('erfStartParams', [2, -1], @isnumeric)
p.addParamValue('plotInconsistentFits', false, @islogical)
p.addParamValue('plotAllFits', false, @islogical)
p.addParamValue('compareToCumGaussFitter', false, @islogical)
p.addParamValue('useRobustFitting', true, @islogical)
p.addParamValue('checkForAnalysisGaps', true, @islogical)

p.parse(elecResp, bootstrapReps, varargin{:})

recalcAll = p.Results.recalcAll;
%keepLogBased = p.Results.keepLogBased;
erfStartParams = p.Results.erfStartParams;
plotInconsistentFits = p.Results.plotInconsistentFits;
plotAllFits = p.Results.plotAllFits;
compareToCumGaussFitter = p.Results.compareToCumGaussFitter;
useRobustFitting = p.Results.useRobustFitting;
checkForGaps = p.Results.checkForAnalysisGaps;

nMovies = length(elecResp.stimInfo.movieNos);

%% checks to make sure that elecResp.analysis.stimAmps reflects the amplitude on the primary electrode
% (due to a change in how stimAmps was calculated in older vs. current versions of code)

if isfield(elecResp.stimInfo, 'pElec')
    trueStimAmp = getStimAmps(elecResp.names.data_path, elecResp.stimInfo.patternNo,...
        elecResp.stimInfo.movieNos(1));
    trueStimAmp = trueStimAmp(elecResp.stimInfo.electrodes == elecResp.stimInfo.pElec);
    if elecResp.stimInfo.stimAmps(1) ~= trueStimAmp
        disp('fixing stimAmps...')
        for j = 1:nMovies
            tempStimAmps = getStimAmps(elecResp.names.data_path, elecResp.stimInfo.patternNo,...
                elecResp.stimInfo.movieNos(j));
            elecResp.stimInfo.stimAmps(j) = tempStimAmps(elecResp.stimInfo.electrodes == elecResp.stimInfo.pElec);
        end
        elecResp.analysis.erfCurrent = 0;
    end
else
    disp('No primary electrode has been specified (elecResp.stimInfo.pElec doesn''t exist).')
end

%% checks to make sure values in successRates and nPulses are correct and finds minimum success rate from
% analyzed movies
minSuccessRate = 1;

for i = 1:nMovies
    mNum = elecResp.stimInfo.movieNos(i);
    if elecResp.stimInfo.nPulses(i) == 0 || recalcAll %value hasn't been set yet
        dataTraces=NS_ReadPreprocessedData(elecResp.names.data_path, '', 0, elecResp.stimInfo.patternNo, mNum, 99999);
        elecResp.stimInfo.nPulses(i) = size(dataTraces, 1);
        elecResp.analysis.erfCurrent = 0;
    elseif ~isempty(elecResp.analysis.type{i}) && elecResp.stimInfo.nPulses(i) ~= length(elecResp.analysis.latencies{i})
        warndlg('number of pulses listed in elecResp.stimInfo.nPulses(i) ~= length of elecResp.analysis.latencies{i}--correcting')
        dataTraces=NS_ReadPreprocessedData(elecResp.names.data_path, '', 0, elecResp.stimInfo.patternNo, mNum, 99999);
        elecResp.stimInfo.nPulses(i) = length(elecResp.analysis.latencies{i});
        if elecResp.stimInfo.nPulses(i) ~= size(dataTraces, 1);
            errordlg('number of pulses in preprocessed data does not match length of elecResp.analysis.latencies{i}')
        end
    end
    if ~isempty(elecResp.analysis.type{i}) %movie has been analyzed
        successRate = sum(elecResp.analysis.latencies{i}~=0)/length(elecResp.analysis.latencies{i});
        if elecResp.analysis.successRates(i) ~= successRate
            disp(['warning: success rate for movie ' num2str(mNum) ' had to be corrected'])
            elecResp.analysis.successRates(i) = successRate;
            elecResp.analysis.erfCurrent = 0;
        end
        minSuccessRate = min(minSuccessRate, elecResp.analysis.successRates(i));
    end
end


%% checks for stimulus amplitudes that haven't been analyzed, other than at
% the beginning or end of the curve

if checkForGaps
    prevNotEmpty = 0;
    prevNotEmptyThenEmpty = 0;
    for j = 1:nMovies
        if isempty(elecResp.analysis.type{j})
            if prevNotEmpty
                prevNotEmptyThenEmpty = 1;
            end
        else
            prevNotEmpty = 1;
            if prevNotEmptyThenEmpty
                disp(['warning: elecResp for ' elecResp.names.rrs_short_name...
                    ' has internal movie numbers that haven''t been analyzed'])
                break
            end
        end
    end
end

%% checks if all analyzed movies are 'locked'
for j=1:nMovies
    if ~isempty(elecResp.analysis.type{j}) && ~elecResp.analysis.finalized(j)
        %try
        %    cprintf([1,0,1], '%s\n', ['warning: not all analyzed movies in ' elecResp.names.rrs_short_name ' have been finalized'])
        %catch
            disp(['warning: not all analyzed movies in ' elecResp.names.rrs_short_name ' have been finalized'])
        %end
        break
    end
end


%% checks if full curve has been analyzed (response probability reaches at least 95%)

% if max(elecResp.analysis.successRates) < 0.95
%     disp(['warning: cell doesn''t reach 0.95 response probability in the analyzed portion of ' elecResp.names.rrs_short_name])
% end

%% recalculates erf if not current (and standard deviation)

data = zeros(3, nMovies);
data(1,:) = abs(elecResp.stimInfo.stimAmps);
data(2,:) = elecResp.analysis.successRates;
data(3,:) = elecResp.stimInfo.nPulses;

for j = length(elecResp.stimInfo.stimAmps): -1: 1
    if isempty(elecResp.analysis.type{j})
        data(:,j) = [];
    end
end

%removes outdated field: gaussCdfParams
if isfield(elecResp.analysis, 'gaussCdfParams')
    elecResp.analysis = rmfield(elecResp.analysis, 'gaussCdfParams');
end

if max(elecResp.analysis.successRates) < 0.2 || minSuccessRate > 0.8
    disp(['warning: cell doesn''t reach 0.2 or starts at higher than 0.8 response probability in the analyzed portion of ' elecResp.names.rrs_short_name ': threshold set to infinity'])
    elecResp.analysis.erfParams = [0 0];
    elecResp.analysis.erfErr = 0;
    if max(elecResp.analysis.successRates) < 0.2
        elecResp.analysis.threshold = inf;
    else %minSuccessRate > 0.8
        elecResp.analysis.threshold = -inf;
    end
    elecResp.analysis.threshStd = 0;
    elecResp.analysis.details.bootstrapReps = 0;

else
    if ~elecResp.analysis.erfCurrent || recalcAll || ~isfield(elecResp.analysis, 'threshStd')
        %recalculate erf fit
        [elecResp.analysis.erfParams , ~, elecResp.analysis.erfErr] = erfFitter(data, erfStartParams(1), erfStartParams(2), 'makePlot', plotAllFits, 'useRobustFitting', useRobustFitting);
        elecResp.analysis.threshold = -elecResp.analysis.erfParams(2)/elecResp.analysis.erfParams(1);
    end
    
    if compareToCumGaussFitter %doesn't store gauss cdf parameters, just checks if fitting is consistent with erf
        %same curve as erf but parameters repackaged
        [gaussCdfParams gaussCdfProjection] = gaussCdfFitter(data, 'makePlot', plotAllFits);
        
        %check if equivalent parameters for erf and cumulative guassian fits are within 0.1% error of
        %eachother
        if abs(gaussCdfParams(1) - elecResp.analysis.threshold) > 0.001*gaussCdfParams(1) ||...
                abs(gaussCdfParams(2) - (-sqrt(2))*elecResp.analysis.erfParams(2)) > 0.001*gaussCdfParams(2)
            disp('warning: parameter values for erf fit and cumulative gaussian fit are not equivalent-- resetting gauss cdf params to match erf params')
            
            if plotInconsistentFits
                %plot erf fit result
                [elecResp.analysis.erfParams , ~, elecResp.analysis.erfErr] = erfFitter(data, erfStartParams(1), erfStartParams(2), 'makePlot', 1, 'useRobustFitting', useRobustFitting);
                
                %plot cumulative gauss fit result on top
                hold on
                plot(gaussCdfProjection(1,:), gaussCdfProjection(2,:),'b-');
                text(0.1, 0.9, ['pattern ' num2str(elecResp.stimInfo.patternNo)], 'units', 'normalized')
                hold off
                drawnow
            end
        end      
    end
    
    
    % calculates standard deviation of threshold if not already done
    % (linear-based erf)
    if ~isfield(elecResp.analysis, 'threshStd') || ~elecResp.analysis.erfCurrent || ~isfield(elecResp.analysis.details, 'bootstrapReps') ||...
            bootstrapReps ~= elecResp.analysis.details.bootstrapReps || recalcAll

        disp(['calculating threshold standard deviation (linear scale) for ' elecResp.names.rrs_short_name])
        stdev = bootstrapThresh(elecResp, bootstrapReps);
        elecResp.analysis.threshStd = stdev;
        elecResp.analysis.details.bootstrapReps = bootstrapReps;
    end
end

%% removes log-based erf analysis if present (obsolete)

if isfield(elecResp.analysis, 'logErfParams')
    elecResp.analysis = rmfield(elecResp.analysis, 'logErfParams');
end
if isfield(elecResp.analysis, 'logErfErr')
    elecResp.analysis = rmfield(elecResp.analysis, 'logErfErr');
end
if isfield(elecResp.analysis, 'logThresh')
    elecResp.analysis = rmfield(elecResp.analysis, 'logThresh');
end
if isfield(elecResp.analysis, 'logThreshStd')
    elecResp.analysis = rmfield(elecResp.analysis, 'logThreshStd');
end
if isfield(elecResp.analysis.details, 'logBootsrapReps')
    elecResp.analysis.details = rmfield(elecResp.analysis.details, 'logBootstrapReps');
end

%%
elecResp.analysis.erfCurrent = 1;