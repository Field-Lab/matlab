function [tp lambdas excludedPairsBin mValidBound] = linFitSimul(thresholds, details, varargin)

%
% fits lines to paired-stim plots simultaneously
% i.e. fits entire linear model:
%     P = f(I_0 + lam_1*I_1 + lam_2*I_2 + ....)
%
% since measurements are thresholds, what's actually being fit is:
% I_0 = t_p - lam_1*I_1 - lam_2*I_2 - ...
%   where t_p and lam_i are being fit, and I_0 are threshold values
%
%
% outsideValReg: binary vector flagging which of the patterns should not be
% used in fitting


p = inputParser;

p.addRequired('thresholds', @isnumeric)
p.addRequired('details', @isstruct)

p.addParamValue('outsideValReg', [], @islogical)
p.addParamValue('fitWithinValidRegions', true, @islogical)

p.parse(thresholds, details, varargin{:})

invalBin = p.Results.outsideValReg;
fitWithinValidRegions = p.Results.fitWithinValidRegions;


nPatterns = length(details.actualRelAmps);
meanRelAmps = cellfun(@(x) mean(x,2), details.actualRelAmps, 'uniformOutput', false);
meanRelAmps = cell2mat(meanRelAmps')';

if isempty(invalBin)
    invalBin = false(1,nPatterns);
end

I_0 = zeros(0,1); %measured thresholds (in terms of primary current amplitude), observation x 1
I_i = zeros(0,6); %secondary electrode current amplitudes, observation x electrode

sAloneThresh = zeros(6,2); %electrode x [cathodal anodal]

excludeBin = false(0,1); %keeps track of which patterns are to be specifically discluded from fitting
toPlotPInd = zeros(0,1); %keeps track of original pattern indeces

%% extract data
for ii = 1:nPatterns
    % primary alone cathodal or electrode pairs
    %if (sum(any(details.sAmps{ii},2) ~= 0) == 1 && any(details.pAmps{ii})) &&...
    %        ~isinf(thresholds(ii)) || (details.pAloneBin(ii) && details.pAmps{ii}(1) < 0)
    if (sum(meanRelAmps(ii,:)~=0) == 1 && any(details.pAmps{ii})) ...
            || (details.pAloneBin(ii) && details.pAmps{ii}(1) < 0)
        
        badSAmps = false;
        if ~details.pAloneBin(ii) %pair
            
            sInd = find(meanRelAmps(ii,:)~=0);
            
            relAmp = meanRelAmps(ii,sInd);
            
            %check if relative amplitude is small enough to have been rounded
            %down to zero at some point in stimulus
            if ~details.pAloneBin(ii) && any(details.sAmps{ii}(sInd,:)==0)
                disp(['pattern index ' num2str(ii) ' was not included in fit because secondary current amplitude sometimes rounded down to 0'])
                badSAmps = true;
            end
            
            I_i = [I_i; zeros(1,6)]; %#ok<AGROW>
            I_i(end, sInd) = relAmp*thresholds(ii);
            
        else
            I_i = [I_i; zeros(1,6)]; %#ok<AGROW>
            
            if exist('pAloneThresh', 'var')
                error('found more than one primary-alone cathodal pattern')
            end
            pAloneThresh = thresholds(ii);
        end
        
        I_0 = [I_0; thresholds(ii)]; %#ok<AGROW>

        excludeBin = [excludeBin; (details.excludeFromFits(ii) || invalBin(ii)) || badSAmps || isinf(thresholds(ii))]; %#ok<AGROW>
        toPlotPInd = [toPlotPInd ii]; %#ok<AGROW>
        
        if isinf(thresholds(ii))
            disp(['pattern index ' num2str(ii) ' was not included in fit because it had an infinite threshold!!'])
        end

                
    % secondary electrode alone thresholds
    elseif details.sAloneBin(ii)
        sInd = find(details.sAmps{ii}(:,1)~=0);
        if details.sAmps{ii}(sInd,1) < 0 %corresponds to secondary-alone cathodal
            sAloneThresh(sInd, 1) = thresholds(ii);
        elseif details.sAmps{ii}(sInd,1) > 0 %corresponds to secondary-alone anodal
            sAloneThresh(sInd, 2) = thresholds(ii);
        end
    end
end


%% do the fitting!

if fitWithinValidRegions %mark values outside limits set by secondary-alone thresholds
    % start by excluding points that fall in "invalid" region for lambda = 0
    % (since line has not been fit)
    
    toFitBin = ~excludeBin;
    
    for ii = 1:6        
        if ~isinf(sAloneThresh(ii,1))
            mLim = pAloneThresh/sAloneThresh(ii,1); %slope of (estimated) limit region edge
            toFitBin(I_i(:,ii)>0 & abs(I_0./I_i(:,ii))<abs(mLim)) = false;
        end
        if ~isinf(sAloneThresh(ii,2))
            mLim = -pAloneThresh/sAloneThresh(2); %slope of (estimated) limit region edge
            toFitBin(I_i(:,ii)<0 & abs(I_0./I_i(:,ii))<abs(mLim)) = false;
        end
    end
    
    % now iteratively perform fit and update toFitBin based on fit values until
    % convergence
    converged = false;
    iCount = 0;    
    
    while ~converged && iCount <= 10
        iCount = iCount+1;

        %set up matrices for least-squares regression
        A = I_i(toFitBin,:);
        A = [A ones(size(A,1), 1)]; %last column is ones for y-intercept

        b = I_0(toFitBin);

        x = A\b;

        %lambdas = -x(1:6); %lambda is the negative of the linear coefficient
        mTmp = x(1:6);
        tpTmp = x(7);

        
        % re-determine which points should be excluded from fitting
        toFitBinPrev = toFitBin;
        toFitBin = ~excludeBin;
        
        mLim = zeros(2,6);
        for ii = 1:6
            if ~isinf(sAloneThresh(ii,1))
                mLim(1,ii) = (mTmp(ii)*sAloneThresh(ii,1)+tpTmp)/sAloneThresh(ii,1); %slope of limit region edge
                toFitBin(I_i(:,ii)>0 & abs(I_0./I_i(:,ii))<abs(mLim(1,ii))) = false;
            end
            if ~isinf(sAloneThresh(ii,2))
                mLim(2,ii) = (mTmp(ii)*sAloneThresh(ii,2)-tpTmp)/sAloneThresh(ii,2); %slope of limit region edge
                toFitBin(I_i(:,ii)<0 & abs(I_0./I_i(:,ii))<abs(mLim(2,ii))) = false;
            end
        end
        
        % check if excluded points are the same in the previous iteration
        converged = isequal(toFitBin, toFitBinPrev);
        
%         %for checking function
%         figure
%         for ii = 1:6
%             subplot(3,2,ii); hold on
%             plot(I_i(I_i(:,ii)~=0 & toFitBin, ii), I_0(I_i(:,ii)~=0 & toFitBin), 'k.')
%             plot(I_i(I_i(:,ii)~=0 & ~toFitBin, ii), I_0(I_i(:,ii)~=0 & ~toFitBin), 'kx')
%             xRange = [min(I_i(:,ii)) max(I_i(:,ii))];
%             plot(xRange, xRange*mTmp(ii)+tpTmp, 'k-')
%             plot([0 xRange(2)], [0 mLim(1,ii)*xRange(2)], 'r-')
%             plot([0 xRange(1)], [0 mLim(2,ii)*xRange(1)], 'r-')
%             axis equal
%             title(num2str(ii))
%         end
%         keyboard
    end
    if ~converged
        error('didn''t converge within 10 iterations!')
    end
end

%%


lambdas = -mTmp; %lambdas are the negatives of the linear coefficients
tp = tpTmp;

excludedPairsBin = false(nPatterns,1);
excludedPairsBin(toPlotPInd(~toFitBin)) = true;

mValidBound = mLim;


