function [tp lambdas excludedPairsBin mBound model] = linFitSimulMinError(thresholds, details, varargin)
    
% fits lines to paired-stim plots simultaneously, including regions where
% secondary electrode behaves as the primary driving electrode
% i.e. fits entire linear model:
%     P = f(I_0 + lam_1*I_1 + lam_2*I_2 + ....)
%
% since measurements are thresholds, what's actually being fit is:
% I_0 = t_p - lam_1*I_1 - lam_2*I_2 - ...
%   where t_p and lam_i are being fit, and I_0 are threshold values



p = inputParser;

p.addRequired('thresholds', @isnumeric)
p.addRequired('details', @isstruct)

%threshold mean squared error (squared distance from point to fit line)
%before additional "regions" are added to account for nonlinearities
p.addParamValue('linErrThresh', [], @isnumeric) 
p.addParamValue('plotFits', false, @islogical)

p.parse(thresholds, details, varargin{:})

linErrThresh = p.Results.linErrThresh;
plotFits     = p.Results.plotFits;


nPatterns = length(details.actualRelAmps);
meanRelAmps = cellfun(@(x) mean(x,2), details.actualRelAmps, 'uniformOutput', false);
meanRelAmps = cell2mat(meanRelAmps')';

excludedPairsBin = false(nPatterns,1); %flags any pair patterns that are excluded from fitting (in any region)

% %initialize structures for storing data and model values
for ii = 1:6
    %stores data points belonging to each of the 6 sElecs  (p-alone thresh
    %will be included in all)
    data(ii).sVal = []; %secondary electrode amplitude at threshold
    data(ii).pVal = []; %primary electrode amplitude at threshold
    data(ii).pAloneThresh = 0;
    data(ii).sAloneThresh = [inf inf]; %positive, negative secondary-alone
%     
%     %model intercepts
%     model(ii).tp = 0; %intercept of model through primary electrode axis
%     model(ii).tsPos = inf; %intercept of model through positive secondary electrode axis (inf signifies no measured sec-alone thresh)
%     model(ii).tsNeg = inf; %intercept of model through positive secondary electrode axis (inf signifies no measured sec-alone thresh)
% 
%     %lambda values: 1 for each region (primary, secondary positive,
%     %secondary negative
%     model(ii).pLam = 0;
%     model(ii).sLamPos = inf;
%     model(ii).sLamNeg = inf;
end


%% extract data
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
                    disp(['pattern index ' num2str(ii) ' was not included in linear model fit because it had an infinite threshold!!'])
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
                
                if details.threshMins(ii) && ~isinf(thresholds(ii))
                    disp('secondary-alone threshold exists for cell that doesn''t reach minimum response probabality')
                    keyboard
                end
            end
            
            if ~isempty(pAmp) && ~isinf(sAmp) && ~excludedPairsBin(ii)
                data(jj).pVal = [data(jj).pVal; pAmp];
                data(jj).sVal = [data(jj).sVal; sAmp];
            end      
        end
    end
end


%% perform fit
%options = optimset('MaxFunEvals', 50000, 'MaxIter', 50000, 'FunValCheck', 'on', 'tolx', 1e-60);
options = optimset('MaxFunEvals', 50000, 'MaxIter', 50000, 'FunValCheck', 'on', 'tolx', 1e-60);

% run an initial fit separately for each secondary electrode that may be
% split into multiple regions, to see if multiple regions are required or
% if it would just be overfitting

pLamStart = zeros(1,6);
if ~isempty(linErrThresh)
    splitIntoRegions = false(1,6);
    for ii = 1:6
        if ~(isinf(data(ii).sAloneThresh(1)) && isinf(data(ii).sAloneThresh(2))) %only bother if there's at least one measured s-alone threshold
            paramsStart(1,1) = 0;
            labels(1).sInd = 1;
            labels(1).fieldName = 'pLam';
            
            paramsStart(2,1) = data(ii).pAloneThresh;
            labels(2).sInd = 1;
            labels(2).fieldName = 'tp';
            
            %remove sAloneThresh from data
            tmpData = data(ii);
            sAloneInd = data(ii).pVal == 0;
            tmpData.pVal(sAloneInd) = [];
            tmpData.sVal(sAloneInd) = [];
            clear noP
            
            [singRegParams, error] = fminsearch(@(paramsFit)linFitSimulError(paramsFit, labels, tmpData, []), paramsStart, options);
            
            if strcmp(labels(1).fieldName, 'pLam')
                pLamStart(ii) = singRegParams(1,1);
            else
                error('unexpected paramsStart structure')
            end
            
            normErr = error/length(data(ii).pVal); %normalize by number of data points
            normErr = normErr/(range(data(ii).sVal)^2 + range(data(ii).pVal)^2); %normalize by range of secondary electrode values in data
            
            disp(['secondary electrode ' num2str(ii) ' nonlinearity index = ' num2str(normErr)])
            if normErr > linErrThresh
                splitIntoRegions(ii) = true;
                disp(['secondary electrode ' num2str(ii) ' qualifies for splitting into regions (exceeds nonlinearity index threshold)'])
            end
        end
    end
    
    
    %for secondary pair data that satisfies linearity condition, remove
    %secondary-alone threshold values from data so that they aren't used in
    %single-region fitting
    for ii = 1:6
        if ~splitIntoRegions(ii)
            sAloneInd = data(ii).pVal == 0;
            data(ii).pVal(sAloneInd) = [];
            data(ii).sVal(sAloneInd) = [];
        end
    end
end


%% 
% %format parameters as a vector of non-infinite values and a struct array
% %describing what each parameter value signifies
% 
% paramInd = 0;
% %initialize slopes
% for ii = 1:6
%     [~, sortOrd] = sort(data(ii).sVal./data(ii).pVal);
%     pAloneInOrd = find(data(ii).sVal(sortOrd)==0);
%     regPtsX = data(ii).sVal(sortOrd([pAloneInOrd-1 pAloneInOrd+1]));
%     regPtsY = data(ii).pVal(sortOrd([pAloneInOrd-1 pAloneInOrd+1]));
%     slopeStart = -diff(regPtsY)/diff(regPtsX);
%     
%     paramInd = paramInd+1;
%     paramsStart(paramInd,1) = slopeStart; %initial slope for primary cathodal region
%     labels(paramInd).sInd = ii;
%     labels(paramInd).fieldName = 'pLam';
%     
%     %only initialize slope for secondary regions if secondary-alone
%     %thresholds exist
%     if ~isinf(data(ii).sAloneThresh(1)) && splitIntoRegions(ii)
%         
%         % determine a reasonable starting lambda for region:        
%         % pick out two points most likely to be in region (most positive
%         % secondary thresholds), calculate slope of line between them,
%         % convert to lambda value for region
%         regPtsX = data(ii).sVal(sortOrd(end-1:end));
%         regPtsY = data(ii).pVal(sortOrd(end-1:end));
%         slopeStart = -diff(regPtsX)/diff(regPtsY);
%         
%         paramInd = paramInd+1;
%         paramsStart(paramInd,1) = slopeStart; %initial lambda for secondary cathodal region
%         labels(paramInd).sInd = ii;
%         labels(paramInd).fieldName = 'sLamPos';
%     end
%     if ~isinf(data(ii).sAloneThresh(2)) && splitIntoRegions(ii)
%         
%         % determine a reasonable starting slope for region:
%         % pick out two points most likely to be in region (most negative
%         % secondary thresholds) and calculate slope of line between them
%         regPtsX = data(ii).sVal(sortOrd(1:2));
%         regPtsY = data(ii).pVal(sortOrd(1:2));
%         slopeStart = -diff(regPtsX)/diff(regPtsY);
%         
%         paramInd = paramInd+1;
%         paramsStart(paramInd,1) = slopeStart; %initial slope for secondary anodal region
%         labels(paramInd).sInd = ii;
%         labels(paramInd).fieldName = 'sLamNeg';
%     end
% end
% 
% %store intercepts separately because they will not be free parameters in
% %initial fitting
% intInd = 0;
% for ii = 1:6
%     %only add one tp value because it is yoked for all secondary electrodes
%     if ii == 1
%         paramInd = paramInd+1;
%         intInd = intInd+1;
%         intercepts(intInd,1) = data(ii).pAloneThresh;
%         labels(paramInd).sInd = ii;
%         labels(paramInd).fieldName = 'tp';
%     end
%     %model(ii).tp = data(ii).pAloneThresh;
%     
%     if ~isinf(data(ii).sAloneThresh(1)) && splitIntoRegions(ii)
%         paramInd = paramInd+1;
%         intInd = intInd+1;
%         intercepts(intInd,1) = data(ii).sAloneThresh(1);
%         labels(paramInd).sInd = ii;
%         labels(paramInd).fieldName = 'tsPos';
%         %model(ii).tsPos = data(ii).sAloneThresh(1);
%     end
%     if ~isinf(data(ii).sAloneThresh(2)) && splitIntoRegions(ii)
%         paramInd = paramInd+1;
%         intInd = intInd+1;
%         intercepts(intInd,1) = data(ii).sAloneThresh(2);
%         labels(paramInd).sInd = ii;
%         labels(paramInd).fieldName = 'tsNeg';
%         %model(ii).tsNeg = data(ii).sAloneThresh(2);
%     end
% end


%% fitting!

%first, fit each electrode pair separately (with intercept locked)
singPairLamFits = cell(1,6);
interceptsAll = cell(1,6);
labelsAll = cell(1,6);
for ii = 1:6
    paramInd = 1;
    tmpData = data(ii);
    
    paramsStart = [];
    
    %initial lambda estimate
     [~, sortOrd] = sort(tmpData.sVal./tmpData.pVal);
%     pAloneInOrd = find(tmpData.sVal(sortOrd)==0);
%     regPtsX = tmpData.sVal(sortOrd([pAloneInOrd-1 pAloneInOrd+1]));
%     regPtsY = tmpData.pVal(sortOrd([pAloneInOrd-1 pAloneInOrd+1]));
%     slopeStart = -diff(regPtsY)/diff(regPtsX);
%     
%     paramsStart(paramInd,1) = slopeStart;
    paramsStart(paramInd,1) = 0;
    labelsAll{ii}(paramInd).sInd = 1;
    labelsAll{ii}(paramInd).fieldName = 'pLam';
    
    if ~isinf(tmpData.sAloneThresh(1)) && splitIntoRegions(ii)
        
        % determine a reasonable starting lambda for region:
        % pick out two points most likely to be in region (most positive
        % secondary thresholds), calculate slope of line between them,
        % convert to lambda value for region
        regPtsX = tmpData.sVal(sortOrd(end-1:end));
        regPtsY = tmpData.pVal(sortOrd(end-1:end));
        slopeStart = -diff(regPtsX)/diff(regPtsY);
        
        paramInd = paramInd+1;
        paramsStart(paramInd,1) = slopeStart;
        labelsAll{ii}(paramInd).sInd = 1;
        labelsAll{ii}(paramInd).fieldName = 'sLamPos';
    end
    if ~isinf(tmpData.sAloneThresh(2)) && splitIntoRegions(ii)
        
        % determine a reasonable starting slope for region:
        % pick out two points most likely to be in region (most negative
        % secondary thresholds) and calculate slope of line between them
        regPtsX = tmpData.sVal(sortOrd(1:2));
        regPtsY = tmpData.pVal(sortOrd(1:2));
        slopeStart = -diff(regPtsX)/diff(regPtsY);
        
        paramInd = paramInd+1;
        paramsStart(paramInd,1) = slopeStart;
        labelsAll{ii}(paramInd).sInd = 1;
        labelsAll{ii}(paramInd).fieldName = 'sLamNeg';
    end
    
    % intercepts (constant)
    intInd = 1;
    
    paramInd = paramInd+1;
    interceptsAll{ii}(intInd,1) = tmpData.pAloneThresh;
    labelsAll{ii}(paramInd).sInd = 1;
    labelsAll{ii}(paramInd).fieldName = 'tp';
    
    if ~isinf(tmpData.sAloneThresh(1)) && splitIntoRegions(ii)
        paramInd = paramInd+1;
        intInd = intInd+1;
        interceptsAll{ii}(intInd,1) = tmpData.sAloneThresh(1);
        labelsAll{ii}(paramInd).sInd = 1;
        labelsAll{ii}(paramInd).fieldName = 'tsPos';
    end
    if ~isinf(tmpData.sAloneThresh(2)) && splitIntoRegions(ii)
        paramInd = paramInd+1;
        intInd = intInd+1;
        interceptsAll{ii}(intInd,1) = tmpData.sAloneThresh(2);
        labelsAll{ii}(paramInd).sInd = 1;
        labelsAll{ii}(paramInd).fieldName = 'tsNeg';
    end
    
%     if ii == 6
%         figure
%         singPairLamFits{ii} = fminsearch(@(lambdas)linFitSimulError([lambdas; interceptsAll{ii}], labelsAll{ii}, tmpData, 1, 'plotAxes', {gca}), paramsStart, options);
%     else
        singPairLamFits{ii} = fminsearch(@(lambdas)linFitSimulError([lambdas; interceptsAll{ii}], labelsAll{ii}, tmpData, []), paramsStart, options);
%    end
end

%combine single pair fit parameters
paramsStart = [];
labels = [];
for ii = 1:6
    paramsStart = [paramsStart; singPairLamFits{ii}; interceptsAll{ii}];
    
    for jj = 1:length(labelsAll{ii})
        labelsAll{ii}(jj).sInd = ii;
    end
    labels = [labels labelsAll{ii}];
end


% refit with intercepts as free parameters

paramsFit = fminsearch(@(params)linFitSimulError(params, labels, data, []), paramsStart, options);

% plot current state of model and get data labels (which stimulus space
% region each datapoint is assigned to)
% note: important to check that points are assigned to regions as expected
% from drawing straight boundaries between regions (going through origin)

if plotFits || any(splitIntoRegions)
    %[~, dataLabels] = linFitSimulError(paramsStart, labels, data, 1:6);
    figure
    plotAxes{1} = axes('position', [0.1 0.1 0.4 0.25]);
    plotAxes{2} = axes('position', [0.55 0.4 0.4 0.25]);
    plotAxes{3} = axes('position', [0.1 0.7 0.4 0.25]);
    plotAxes{4} = axes('position', [0.55 0.1 0.4 0.25]);
    plotAxes{5} = axes('position', [0.1 0.4 0.4 0.25]);
    plotAxes{6} = axes('position', [0.55 0.7 0.4 0.25]);
    [~, dataLabels] = linFitSimulError(paramsFit, labels, data, 1:6, 'plotAxes', plotAxes);
    
    for ii = 1:6
        if splitIntoRegions(ii)
            set(plotAxes{ii}, 'color', [0.9 0.9 0.9])
        end
    end
    
else
    [~, dataLabels] = linFitSimulError(paramsFit, labels, data, []);
end


%% unpack final parameters

%initialize structures for storing data and model values
for ii = 1:6
    %model intercepts
    model(ii).tp = 0; %#ok<*AGROW> %intercept of model through primary electrode axis
    model(ii).tsPos = inf; %intercept of model through positive secondary electrode axis (inf signifies no measured sec-alone thresh)
    model(ii).tsNeg = inf; %intercept of model through positive secondary electrode axis (inf signifies no measured sec-alone thresh)

    %lambda values: 1 for each region (primary, secondary positive,
    %secondary negative
    model(ii).pLam = 0;
    model(ii).sLamPos = inf;
    model(ii).sLamNeg = inf;
end

for ii = 1:length(paramsFit)
    switch labels(ii).fieldName
        case 'tp'
            for jj = 1:length(model)
                model(jj).tp = paramsFit(ii); %same for all secondary electrodes
            end
        case 'tsPos'
            model(labels(ii).sInd).tsPos   = paramsFit(ii);
        case 'tsNeg'
            model(labels(ii).sInd).tsNeg   = paramsFit(ii);
        case 'pLam'
            model(labels(ii).sInd).pLam    = paramsFit(ii);
        case 'sLamPos'
            model(labels(ii).sInd).sLamPos = paramsFit(ii);
        case 'sLamNeg'
            model(labels(ii).sInd).sLamNeg = paramsFit(ii);
        otherwise
            error('unrecognized field name')
    end
end

tp = model(1).tp;
lambdas = [model.pLam];

%% calculate boundaries between regions, defined by slope of line passing
% through origin

%define boundaries between regions of different linear models IFF at least 2
%points are fit in a particular secondary electrode region (and a
%secondary-alone threshold exists)

mBound = zeros(2,6); %slope of lines dividing regions

for ii = 1:6
    %positive primary region: 2 points on line, [sVal; pVal]
    p1 = [0; model(ii).tp];
    p2 = [1; model(ii).tp - model(ii).pLam];
        
    if ~isinf(model(ii).tsPos) && sum(dataLabels(ii).sPos) >= 2
        %find intersection point between lines defined by point pairs
        
        %positive secondary region: 2 points on line [sVal; pVal]
        q1 = [model(ii).tsPos;                     0];
        q2 = [model(ii).tsPos - model(ii).sLamPos; 1];
        
        x = [p1(1) q1(1); p2(1), q2(1)]; %starting points in first row, ending points in second row
        y = [p1(2) q1(2); p2(2), q2(2)];
        
        %algorithm from stackoverflow
        dx = diff(x);  % Take the differences down each column
        dy = diff(y);
        den = dx(1)*dy(2)-dy(1)*dx(2);  % Precompute the denominator
        u = (dx(2)*(y(1)-y(3))-dy(2)*(x(1)-x(3)))/den;
        
        %intersection point
        xi = x(1)+u*dx(1);
        yi = y(1)+u*dy(1);
        
        mBound(1,ii) = yi/xi;
        model(ii).posInt = [xi yi];
    else
        model(ii).posInt = [];
    end
    if ~isinf(model(ii).tsNeg) && sum(dataLabels(ii).sNeg) >= 2
        
        %negative secondary regions: 2 points on line [sVal; pVal]
        q1 = [-model(ii).tsNeg;                     0];
        q2 = [-model(ii).tsNeg - model(ii).sLamNeg; 1];
        
        x = [p1(1) q1(1); p2(1), q2(1)]; %starting points in first row, ending points in second row
        y = [p1(2) q1(2); p2(2), q2(2)];
        
        %algorithm from stackoverflow
        dx = diff(x);  % Take the differences down each column
        dy = diff(y);
        den = dx(1)*dy(2)-dy(1)*dx(2);  % Precompute the denominator
        u = (dx(2)*(y(1)-y(3))-dy(2)*(x(1)-x(3)))/den;
        
        %intersection point
        xi = x(1)+u*dx(1);
        yi = y(1)+u*dy(1);
        
        mBound(2,ii) = yi/xi;
        model(ii).negInt = [xi yi];
    else
        model(ii).negInt = [];
    end

end



            
            
            
            
            
            
            
            
            
            