function [thresholds threshStds details] = extractPatternStimAnalysis(pathToData, patternNos, pElec, neuronID, varargin)
%
%
% arguments:
%   pathToData  should be path to elecResp files, which must be named in
%   the format: 'elecResp_n [neuronID]  _p [pattern number]'
%   patternNos  vector of integers corresponding to pattern numbers in elecResp names
%   pElec       single integer specifying primary electrode
%   neuronID    single integer specifying cell ID (as specified in elecResp file names)
%
% optional parameters:
%   useRobustFitting
%
%
% CHANGES
% 2010-04-07
% changed calculation of estimated primary-alone threshold to reflect better linear fit in terms of
% absolute (vs relative) amplitude on secondary
% 
% 2013-01-07
% automatically sets threshold to inf and erf fit related values to 0 (in addition to calculating a
% threshold minimum) if response probability doesn't reach
% requiredSuccessRate
% *** saved elecResp file is not altered -- only function output is altered
% *** has an effect on figures (formerly plotted thresholds for response
% curves as long as they had probabilities above 0.2)

p = inputParser;

p.addRequired('pathToData', @ischar)
p.addRequired('patternNos', @isnumeric)
p.addRequired('pElec', @isnumeric)
p.addRequired('neuronID', @isnumeric)

p.addParamValue('useRobustFitting', true, @islogical)
p.addParamValue('requiredSuccessRate', 0.4, @isnumeric) %calculates a "minimum possible thresholds" if there are no measured success rates above this value
p.addParamValue('checkForAnalysisGaps', true, @islogical) %whether to display a message when a responds curve has internal current amplitudes that aren't analyzed
p.addParamValue('sElecs', [], @isnumeric) %overrides default of nearest-neighbor electrodes
p.addParamValue('refitAllErf', false, @islogical)
p.addParamValue('excludeFromFits', [], @isnumeric) %flags corresponding patterns
p.addParamValue('plotBSResults', false, @islogical) %if true and if bootstrap analysis is performed, generated distribution of thresholds will be plotted

p.parse(pathToData, patternNos, pElec, neuronID, varargin{:})

useRobustFitting = p.Results.useRobustFitting;
reqSuccRate = p.Results.requiredSuccessRate;
checkForAnalysisGaps = p.Results.checkForAnalysisGaps;
sElecs = p.Results.sElecs;
refitAllErf = p.Results.refitAllErf;
pExclude = p.Results.excludeFromFits;
plotBSResults = p.Results.plotBSResults;

prevDir = pwd;
cd(pathToData)

nPatterns = length(patternNos);

if isempty(sElecs) %not specified as an argument
    allElecs = getCluster(pElec); %first is primary electrode
    sElecs = allElecs(2:end);
    clear allElecs
%else
%    allElecs = [pElec sElecs];
end

%elecs =                  cell(nPatterns, 1); %electrodes used in pattern (phasing out)
%stimElecs =              cell(nPatterns, 1); %logical vectors defining which of allElecs is in current movie has nonzero amplitude pulses for at least 1 movie
%connElecs =              cell(nPatterns, 1); %logical vectors defining which of allElecs has "connect" signal on (but not necessarily a nonzero amplitude pulse)

%baseAmps =               cell(nPatterns, 1); %corresponding electrode amplitudes for first movie
thresholds =             zeros(nPatterns, 1);
threshStds =             zeros(nPatterns, 1);
threshMins =             zeros(nPatterns, 1);
curveWidths =            zeros(nPatterns, 1);

%pAmps =                  zeros(nPatterns, 1);
pAmpsAll =               cell(nPatterns, 1);

%sAmps =                  cell(nPatterns, 1); % each cell is a vector
%corresponding to sElecs
sAmpsAll =               cell(nPatterns, 1);

actualRelAmps =          cell(nPatterns, 1); %actual relative amplitude of secondary electrodes for first movie
%actualRelAmpsAllMov =    cell(nPatterns, 1); %actual relative amplitude of secondary electrodes for all movies (set to infinity if pAmps == 0)
pAloneBin =              false(nPatterns, 1); % true only if no secondary electrodes had nonzero amplitudes or were connected (should only be true for 1 or 2 pattern)
sAloneBin =              false(nPatterns, 1); % true only if nonzero amplitude on a secondary and zero amplitude on primary
sConnected =             cell(nPatterns, 1); % whether "connect" signal was on for any secondary electrodes

excludeFromFits = false(nPatterns, 1);

for ii = 1:nPatterns
    fprintf('%c', '.')
    if mod(ii,100)==0
        fprintf('\n')
    end
    
    fileName = ['elecResp_n' num2str(neuronID) '_p' num2str(patternNos(ii))];
    temp = load(fileName);
    
    if any(pExclude == patternNos(ii))
        excludeFromFits(ii) = true;
    end
    
    disp(num2str(temp.elecResp.analysis.details.bootstrapReps))
    
    elecResp = checkForUnfinishedAnalysis(temp.elecResp, 100, 'recalcAll', false,...
        'useRobustFitting', useRobustFitting, 'checkForAnalysisGaps', checkForAnalysisGaps, 'recalcAll', refitAllErf, 'plotBSResults', plotBSResults);
    
    if elecResp.stimInfo.pElec ~= pElec;
        disp(['primary electrode specified in ' fileName ' does not match primary electrode specified as function argument'])
        keyboard
    end
    
    nMovies = length(elecResp.stimInfo.movieNos);

    %extracts info necessary to plot
    %elecs{ii} = elecResp.stimInfo.electrodes;
    %baseAmps{ii} = getStimAmps(elecResp.names.data_path,
    %elecResp.stimInfo.patternNo, elecResp.stimInfo.movieNos(1));
    
    
    %% extract stimulus information

    if ~isfield(elecResp.stimInfo, 'clusterStimDetails')
        % extract stim info from pattern file
        
        sAmpsAll{ii} = zeros(length(sElecs), nMovies);
        pAmpsAll{ii} = zeros(1, nMovies);
        sConnTmp = false(length(sElecs), nMovies);
        
        for jj = 1:nMovies
            [stimAmpsTmp, stimElecTmp, ~, elecConnTmp] = getStimAmps(elecResp.names.data_path, elecResp.stimInfo.patternNo, elecResp.stimInfo.movieNos(jj));
            
            for kk = 1:length(stimElecTmp)
                if stimElecTmp(kk) == pElec
                    pAmpsAll{ii}(jj) = stimAmpsTmp(kk);
                elseif any(sElecs == stimElecTmp(kk))
                    sAmpsAll{ii}(sElecs == stimElecTmp(kk), jj) = stimAmpsTmp(kk);
                else
                    error(['stimulation occurs on an unexpected electrode in ' fileName])
                end
            end
            for kk = 1:length(elecConnTmp)
                if any(sElecs == elecConnTmp(kk))
                    sConnTmp(sElecs == elecConnTmp(kk), jj) = true;
                elseif elecConnTmp(kk) ~= pElec
                    error(['unexpected electrode in ' fileName ' is ''connected'''])
                end
            end
        end
        
        if any(mod(mean(sConnTmp, 2),1)) %checks to make sure connect signal is the same for all movies
            disp(['connected electrodes vary with movie for ' fileName])
            keyboard
        else
            sConnected{ii} = sConnTmp(:,1); %they're all the same so just store values from first movie
        end
        
        %save secondary pulse amplitudes in elecResp file
        elecResp.stimInfo.clusterStimDetails.sElecs = sElecs;
        elecResp.stimInfo.clusterStimDetails.pAmpsAll = pAmpsAll{ii};
        elecResp.stimInfo.clusterStimDetails.sAmpsAll = sAmpsAll{ii};
        elecResp.stimInfo.clusterStimDetails.sConnected = sConnected{ii};
    else
        pAmpsAll{ii} = elecResp.stimInfo.clusterStimDetails.pAmpsAll;
        sAmpsAll{ii} = elecResp.stimInfo.clusterStimDetails.sAmpsAll;
        sConnected{ii} = elecResp.stimInfo.clusterStimDetails.sConnected;
    end
    
    save(fileName, 'elecResp')
    
    clear fileName
    
    %amplitude on secondary(ies)
    actualRelAmps{ii} = zeros(length(sElecs), nMovies);
    if any(pAmpsAll{ii}) %pulse on primary
        for jj = 1:length(sElecs)
            actualRelAmps{ii}(jj,:) = sAmpsAll{ii}(jj,:)./pAmpsAll{ii};
        end
    elseif any(any(sAmpsAll{ii}))
        sAloneBin(ii) = true;
        for jj = 1:length(sElecs)
            if any(sAmpsAll{ii}(jj,:))
                actualRelAmps{ii}(jj,:) = inf;
            end
        end
    end
    if any(pAmpsAll{ii}) && ~any(any(sAmpsAll{ii})) && ~any(sConnected{ii})
        pAloneBin(ii) = true;
    end
    
    %% extract analysis values
    
    thresholds(ii) = elecResp.analysis.threshold;
    threshStds(ii) = elecResp.analysis.threshStd;
    curveWidths(ii) = -sqrt(2)*elecResp.analysis.erfParams(2);
    
    % determines the minumum possible threshold for incomplete curves
    % (if reqSuccRate not met)
    data = zeros(2, length(elecResp.stimInfo.stimAmps));
    data(1,:) = abs(elecResp.stimInfo.stimAmps);
    data(2,:) = elecResp.analysis.successRates;
    data(3,:) = elecResp.stimInfo.nPulses;
    for jj = length(elecResp.stimInfo.stimAmps): -1: 1
        if isempty(elecResp.analysis.type{jj})
            data(:,jj) = [];
        end
    end
    
    if ~any(data(2,:) >= reqSuccRate)
        minThreshData = data;
        for zz = 1:5
            x = size(minThreshData,2);
            minThreshData(1,x+1) = minThreshData(1,x)*1.1;
            minThreshData(2,x+1) = 1;
            minThreshData(3,x+1) = mean(data(3,:));
        end
        p = erfFitter(minThreshData, 2, -1, 'makePlot', 0);
        threshMins(ii) = -p(2)/p(1);
        
        % didn't meet reqSuccRate, so delete threshold and other erf fit
        % parameters
        if ~isinf(thresholds(ii))
            thresholds(ii) = inf;
            threshStds(ii) = 0;
            curveWidths(ii) = 0;
            disp(['warning: resp. prob. doesn''t reach ' num2str(reqSuccRate) ' in the analyzed portion of ' elecResp.names.rrs_short_name ': threshold set to infinity (outside of elecResp file)'])
            %keyboard
        end
    end
end


%%

details.pAloneBin = pAloneBin;
%details.sAmps = sAmps;
%details.pAmps = pAmps;
%details.actualRelAmps = actualRelAmps;

details.sAmps = sAmpsAll;
details.pAmps = pAmpsAll;
details.actualRelAmps = actualRelAmps;
details.sAloneBin = sAloneBin;
details.curveWidths = curveWidths;
details.sConnected = sConnected;
details.threshMins = threshMins;
details.sElecs = sElecs;
details.excludeFromFits = excludeFromFits;

cd(prevDir)

