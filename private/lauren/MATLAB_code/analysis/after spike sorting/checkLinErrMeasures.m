function checkLinErrMeasures(cellInfo, linErrThresh, varargin)


p = inputParser;

p.addRequired('cellInfo', @isstruct)
p.addRequired('linErrThresh', @isnumeric)


p.addParamValue('plotFits', false, @islogical)

p.parse(cellInfo, linErrThresh, varargin{:})

plotFits = p.Results.plotFits;

%%

%labels assigned by visual inspection
for x = 1 %for code folding
    linFitGood = [1 1 1 1 1 1;
        1 1 1 1 1 1;
        1 1 1 1 1 1;
        1 1 1 1 1 1;
        0 1 0 1 1 0;
        1 1 1 1 1 0;
        1 0 1 1 0 1;
        1 1 1 1 1 1;
        0 0 1 1 0 1;
        1 1 1 1 1 1;
        1 0 1 1 1 1;
        1 0 1 1 0 0;
        1 0 1 1 0 1;
        1 0 1 1 0 0;
        1 1 1 1 1 0];
    
    linFitBad = [0 0 0 0 0 0;
        0 0 0 0 0 0;
        0 0 0 0 0 0;
        0 0 0 0 0 0;
        1 0 1 0 0 0;
        0 0 0 0 0 0;
        0 1 0 0 1 0;
        0 0 0 0 0 0;
        0 0 0 0 0 0;
        0 0 0 0 0 0;
        0 0 0 0 0 0;
        0 0 0 0 0 0;
        0 0 0 0 0 0;
        0 0 0 0 0 0;
        0 0 0 0 0 1];
    
    linFitBorderline = [0 0 0 0 0 0;
        0 0 0 0 0 0;
        0 0 0 0 0 0;
        0 0 0 0 0 0;
        0 0 0 0 0 1;
        0 0 0 0 0 1;
        0 0 0 0 0 0;
        0 0 0 0 0 0;
        1 1 0 0 1 0;
        0 0 0 0 0 0;
        0 1 0 0 0 0;
        0 1 0 0 1 1;
        0 1 0 0 1 0;
        0 1 0 0 1 1;
        0 0 0 0 0 0];
end

singRegFitErr = zeros(15,6);

%%
for kk = 1:15
    thresholds = cellInfo(kk).thresholds;
    details = cellInfo(kk).details;
       
    nPatterns = length(details.actualRelAmps);
    meanRelAmps = cellfun(@(x) mean(x,2), details.actualRelAmps, 'uniformOutput', false);
    meanRelAmps = cell2mat(meanRelAmps')';
    
    excludedPairsBin = false(nPatterns,1); %flags any pair patterns that are excluded from fitting (in any region)
    
    %initialize structure for storing data
    for ii = 1:6
        data(ii).sVal = []; %secondary electrode amplitude at threshold
        data(ii).pVal = []; %primary electrode amplitude at threshold
        data(ii).pAloneThresh = 0;
        data(ii).sAloneThresh = [inf inf]; %positive, negative secondary-alone
    end
    
    % extract data
    for jj = 1:6 %loop through secondary electrodes
        for ii = 1:nPatterns
            pAmp = []; sAmp = [];
            if sum(meanRelAmps(ii,:)~=0) < 2 %not a triplet pattern
                if meanRelAmps(ii,jj) && ~details.sAloneBin(ii) %pair pattern involving jjth sElec
                    badSAmps = false;
                    relAmp = meanRelAmps(ii,jj);
                    
                    pAmp = thresholds(ii);
                    sAmp = thresholds(ii)*relAmp;
                    
                    %exclude a datapoint corresponding to a pair if specified in
                    %details.excludeFromFits, if secondary electrode amplitude
                    %sometimes rounds to zero, or if there was no measured
                    %threshold
                    
                    %check if relative amplitude is small enough to have been rounded
                    %down to zero at some point in stimulus
                    if any(details.sAmps{ii}(jj,:)==0)
                        disp(['pattern index ' num2str(ii) ' was not included in fit because secondary current amplitude sometimes rounded down to 0'])
                        badSAmps = true;
                    end
                    if isinf(thresholds(ii));
                        disp(['pattern index ' num2str(ii) ' of cell index ' num2str(kk) ' was not included in linear model fit because it had an infinite threshold!!'])
                    end
                    if details.excludeFromFits(ii)
                        disp(['pattern index ' num2str(ii) ' was not included in linear model fit because it was marked to be excluded in cell_list_cluster_stim!!'])
                    end
                    excludedPairsBin(ii) = details.excludeFromFits(ii) || badSAmps || isinf(thresholds(ii));
                elseif details.pAloneBin(ii) && details.pAmps{ii}(1) < 0 %primary alone, cathodal
                    pAmp = thresholds(ii);
                    sAmp = 0;
                    data(jj).pAloneThresh = thresholds(ii);
                    
                elseif details.sAloneBin(ii) && details.sAmps{ii}(jj,1) %secondary alone
                    pAmp = 0;
                    sAmp = -1*thresholds(ii)*sign(details.sAmps{ii}(jj,1)); %cathodal threshold = positive sElec axis, anodal = negative
                    if details.sAmps{ii}(jj,1) < 0
                        data(jj).sAloneThresh(1) = thresholds(ii);
                    else
                        data(jj).sAloneThresh(2) = thresholds(ii);
                    end
                end
                
                if ~isempty(pAmp) && ~isinf(sAmp) && ~excludedPairsBin(ii)
                    data(jj).pVal = [data(jj).pVal; pAmp];
                    data(jj).sVal = [data(jj).sVal; sAmp];
                end
            end
        end
    end
    
    
    % perform fit
    options = optimset('MaxFunEvals', 50000, 'MaxIter', 50000, 'FunValCheck', 'on', 'tolx', 1e-20);
    
    
    % run an initial fit separately for each secondary electrode that may be
    % split into multiple regions, to see if multiple regions are required or
    % if it would just be overfitting
    
    normErr1 = zeros(1,6);
    normErr2 = zeros(1,6);
    
    if plotFits
        figure('position', [100 100 800 800])
    end
    for ii = 1:6
        %if ~(isinf(data(ii).sAloneThresh(1)) && isinf(data(ii).sAloneThresh(2))) %only bother if there's at least one measured s-alone threshold
            paramsStart(1,1) = 0;
            labels(1).sInd = 1;
            labels(1).fieldName = 'pLam';
            
            paramsStart(2,1) = data(ii).pAloneThresh;
            labels(2).sInd = 1;
            labels(2).fieldName = 'tp';
            
            %remove sAloneThresh from data
            tmpData = data(ii);
            noP = tmpData.pVal == 0;
            tmpData.pVal(noP) = [];
            tmpData.sVal(noP) = [];
            clear noP
            
            
            [paramsFit error] = fminsearch(@(paramsFit)linFitSimulError(paramsFit, labels, tmpData, []), paramsStart, options);
            
            normErr = error/length(data(ii).pVal); %normalize by number of data points
            %normErr2(ii) = normErr/range(data(ii).sVal);
            normErr2(ii) = normErr/(range(data(ii).sVal)^2 + range(data(ii).pVal)^2); %normalize by range of secondary electrode values in data
            
            if plotFits
                %plot fit
                subplot(3,2,ii)
                linFitSimulError(paramsFit, labels, tmpData, 1, 'plotAxes', {gca});
                title(['cell index ' num2str(kk) ', electrode ' num2str(ii) ' mean normalized error = ' num2str(normErr1(ii)), 10,...
                    ' normalized by squared distance = ' , num2str(normErr2(ii))])
            end
            
            paramsFitAll{kk} = paramsFit;
            
            
        %end
        %keyboard
    end
    
    singRegFitErr(kk,:) = normErr2;
end


%%
% %make histograms of different categories
% singRegFitErr = reshape(singRegFitErr',1,[]);
% linFitGood = logical(reshape(linFitGood',1,[]));
% linFitBad = logical(reshape(linFitBad',1,[]));
% linFitBorderline = logical(reshape(linFitBorderline',1,[]));
% 
% logErr = log(singRegFitErr);
% 
% keyboard
% 
% binEdges = -15:0.2:0;
% if any(logErr > max(binEdges) | logErr < min(binEdges))
%     error('chosen bin edges don''t span all data')
% end
% %binEdges = [-inf binEdges inf];
% 
% hist(logErr, 100)

%%

% figure
% subplot(3,1,1)
% [N, bin] = histc(logErr(linFitGood), binEdges);
% bar(binEdges, N, 'histc')
% 
% subplot(3,1,2)
% [N, bin] = histc(logErr(linFitBorderline), binEdges);
% bar(binEdges, N, 'histc')
% ylabel('number of secondary electrodes')
% 
% subplot(3,1,3)
% [N, bin] = histc(logErr(linFitBad), binEdges);
% bar(binEdges, N, 'histc')
% xlabel('log(normalized error)')

%%
% figure
% [N, bin] = histc(logErr, binEdges);
% %bar(binEdges, N, 'histc')
% %hold on
% %plot(log(0.00025)*[1 1], [0 10])

%% nicer histogram for making figure
logErr = log(reshape(singRegFitErr',1,[]));

binEdges = -15:0.2:0;
if any(logErr > max(binEdges) | logErr < min(binEdges))
    error('chosen bin edges don''t span all data');
end

[N, ~] = histc(logErr, binEdges);

patchVals = [binEdges(1); 0];
for ii = 1:length(binEdges)-1
    patchVals = [patchVals, [binEdges(ii) binEdges(ii+1); N(ii) N(ii)]];
end
patchVals = [patchVals, [binEdges(end); 0]];

figure('position', [200 200 400 200]); hold on
patch(patchVals(1,:), patchVals(2,:), 'r')
plot(log(linErrThresh)*[1 1], [0 10], 'k--')
%plot(linErrThresh*[1 1], [0 10], 'k--')

%chosen examples for figure
%piece-wise linear s elecs
multRegSElecs = [25 27 38 41 90];
%exampleSElecs = [25 56 19 53];
exampleSElecs = [25 56 10 53];

n800SElecs = 13:18;

for ii = 1:length(multRegSElecs)
    plot(logErr(multRegSElecs(ii)), 0, 'g+')
end
for ii = 1:length(exampleSElecs)
    plot(logErr(exampleSElecs(ii)), 0, 'b+')
end
for ii = 1:length(n800SElecs)
    plot(logErr(n800SElecs(ii)), 10, 'm+')
end

%% paired plots for examples

figure('position', [100 100 400 400]*2)
axesPos = {[30 220 160 80]*2, [30 120 160 80]*2, [30 20 160 80]*2, [230 220 160 80]*2};
% neuronInds = [5 10 4 9];
% sElecInds = [1 2 1 5];
% xLimAll = {[-0.5 0.5], [-2.2 2.2], [-4 4],  [-2 1]};
% yLimAll = {[0 0.5],    [0 2.2],    [0 4],   [0 1.5]};

neuronInds = [5 10 2 9];
sElecInds = [1 2 4 5];
xLimAll = {[-0.5 0.5], [-1.8 1.8], [-0.6 0.8],  [-2 1]};
yLimAll = {[0 0.5],    [0 1.8],    [0 0.7],   [0 1.5]};


for jj = 1:4
    detailsCopy = cellInfo(neuronInds(jj)).details;
    
    %make a copy of details with single-region model fitting if necessary
    if any(detailsCopy.mValidBound(:,sElecInds(jj)));
        model = detailsCopy.multiRegModel(:,sElecInds(jj));
        
        model.tp = paramsFitAll{neuronInds(jj)}(2);
        model.tsPos = inf;
        model.tsNeg = inf;
        model.pLam = paramsFitAll{neuronInds(jj)}(1);
        model.sLamPos = inf;
        model.sLamNeg = inf;
        model.posInt = [];
        model.negInt = [];
        
        detailsCopy = cellInfo(neuronInds(jj)).details;
        detailsCopy.multiRegModel(:,sElecInds(jj)) = model;
        detailsCopy.mValidBound(:,sElecInds(jj)) = [0 0];
    end
    
    
    axes('units', 'pixels', 'position', axesPos{jj})
    
    
    plotLinearityDataVsModel(cellInfo(neuronInds(jj)).thresholds, cellInfo(neuronInds(jj)).threshStds, detailsCopy, sElecInds(jj),...
        'showSecAloneThresh', false,...
        'xPlotLim', xLimAll{jj}, 'yPlotLim', yLimAll{jj}, 'plotUnfitPairs', false, 'shadeRegType', 'bootstrapSD');
    
    set(gca, 'xlim', xLimAll{jj}, 'ylim', yLimAll{jj})
    title(['electrode ' num2str(details.sElecs(jj))])
    xlabel('secondary current amplitude')
    ylabel('primary current amplitude')
    
    
end



%%








