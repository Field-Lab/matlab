function cellInfo = plotLambdaVAxonDir(cellInfo, varargin)

% cellInfo should be a struct array with the following fields:
%
%	.pathToAxonDirData = '/snle/lab/Experiments/Array/Analysis/2008-08-26-0/data007-NW/axonDirData.mat';
%   .details.lambdas
%   .details.sElecs %with order corresponding to lambdas


p = inputParser;

p.addRequired('cellInfo', @isstruct)

p.addParamValue('excludeBadElecs', true, @islogical)
p.addParamValue('range90', false, @islogical) %whether or not to reflect angles >90 across the y axis
p.addParamValue('binSize', 30, @isnumeric) %in degrees
p.addParamValue('excludeCellsInd', [], @isnumeric) %index (within cellInfo) of particular cell(s) to exclude from the analysis
p.addParamValue('reflectData', false, @islogical) %whether to duplicate data points onto quadrants 3 and 4, for more symmetric visualization
p.addParamValue('nShuffleReps', 100, @isnumeric)
p.addParamValue('lamLimPlot', 1, @isnumeric) %how scale axes and grid to fit up to this value of lambda

p.parse(cellInfo, varargin{:})

params = p.Results;

if ~isfield(cellInfo(1), 'sAlignedAng')
    
    for ii = 1:length(cellInfo)
        % load measured axon directions
        load(cellInfo(ii).pathToAxonDirData)
        axon_dirs = axonDirData.direction;
        
        
        % calculate mean axon direction
        dir_vect = exp(axon_dirs*(pi/180)*1i);
        mean_vect = mean(dir_vect);
        if abs(mean_vect) < 1e-12; %if length of mean vector is very small, angle is unreliable
            error('couldn''t calculated a mean axon angle')
        end
        
        meanAxAng = atan2(imag(mean_vect), real(mean_vect))*180/pi;
        
        cellInfo(ii).meanAxAng = meanAxAng;
        
        %     %for testing only
        %     figure; hold on
        %     for jj = 1:length(axon_dirs)
        %         plot([0 cosd(axon_dirs(jj))], [0 sind(axon_dirs(jj))], 'k-')
        %     end
        %     plot(cosd(meanAxAng), sind(meanAxAng), 'ro')
        %     axis equal; set(gca, 'xlim', [-1.2 1.2], 'ylim', [-1.2 1.2])
        
        %get angles to sElecs in electrode coordinates used for axon direction
        %measurement
        xCoords = axonDirData.elecCoords.x;
        yCoords = axonDirData.elecCoords.y;
        
        %shift coordinates to be centered in primary electrode
        xCoords = xCoords - xCoords(cellInfo(ii).pElec);
        yCoords = yCoords - yCoords(cellInfo(ii).pElec);
        
        %determine angles of sElecs
        sElecXCoords = xCoords(cellInfo(ii).details.sElecs);
        sElecYCoords = yCoords(cellInfo(ii).details.sElecs);
        
        sAngles = atan2(sElecYCoords, sElecXCoords)*180/pi;
        
        %rotate sElec angles so that the axon is at 0 degrees
        sAligned = sAngles - meanAxAng;
        
        %for testing only
        %     figure
        %     hold on
        %     sColors = lines(6);
        %     for jj = 1:6
        %         plot([0 cosd(sAngles(jj))], [0 sind(sAngles(jj))], '-', 'color', sColors(jj,:), 'linewidth', 2)
        %         plot([0 cosd(sAligned(jj))], [0 sind(sAligned(jj))], '--', 'color', sColors(jj,:), 'linewidth', 2)
        %         text(cosd(sAngles(jj)), sind(sAngles(jj)), num2str(cellInfo(ii).details.sElecs(jj)))
        %     end
        %     plot(cosd(meanAxAng), sind(meanAxAng), 'ko'); axis equal
        %     keyboard
        
        cellInfo(ii).sAlignedAng = sAligned;
    end
end


%% gather lambdas

s_angles_all = [];
lam_all      = [];

for ii = 1:length(cellInfo)
    if ~any(params.excludeCellsInd == ii)
        lam = cellInfo(ii).details.lambdas;
        sAngles = cellInfo(ii).sAlignedAng;
        
        if params.excludeBadElecs
            lam = lam(~ismember(cellInfo(ii).details.sElecs, cellInfo(ii).badSElecs));
            sAngles = sAngles(~ismember(cellInfo(ii).details.sElecs, cellInfo(ii).badSElecs));
        end
        
        s_angles_all = [s_angles_all sAngles]; %#ok<AGROW>
        lam_all      = [lam_all lam]; %#ok<AGROW>
    end
end

%maxlam = max(abs(lam_all));
%lam_all = lam_all/maxlam;

%% plot results

%get all angles between 0 and 180
s_angles_all = rem(s_angles_all, 360); %get everything within -360:360;
s_angles_all(s_angles_all<0) = s_angles_all(s_angles_all<0)+360; %get everything within 0:360
s_angles_all(s_angles_all>180) = 360 - s_angles_all(s_angles_all>180); %mirror across x-axis

if params.range90
    s_angles_all(s_angles_all>90) = 180 - s_angles_all(s_angles_all>90); %mirror across y-axis
end



%% shuffled angle analysis
    
if params.range90
    binLims = 0:params.binSize:90;
else
    binLims = 0:params.binSize:180;
end
nBins = length(binLims)-1;

binMeansAllShuffled = inf*ones(params.nShuffleReps, nBins);
binCountAllShuffled = zeros(params.nShuffleReps, nBins);

%sum of squared residuals measure

residAllShuff = zeros(params.nShuffleReps, 1);
%entropy measures
entAllShuff = zeros(params.nShuffleReps, 1);

%peak value measure
peakAllShuff = zeros(params.nShuffleReps, 1);

lamMean = mean(lam_all);

tic
for kk = 1:params.nShuffleReps
    s_angles_shuffled = [];
    lam_all_shuff      = [];
    
    for ii = 1:length(cellInfo)
        randShift = 360*rand(1); %add a random angle between 0 and 360 to each cell
        
        if ~any(params.excludeCellsInd == ii)
            lam = cellInfo(ii).details.lambdas;
            sAngles = cellInfo(ii).sAlignedAng + randShift;
            
            if params.excludeBadElecs
                lam = lam(~ismember(cellInfo(ii).details.sElecs, cellInfo(ii).badSElecs));
                sAngles = sAngles(~ismember(cellInfo(ii).details.sElecs, cellInfo(ii).badSElecs));
            end
            
            
%             %get all angles between 0 and 180
%             tmp = rem(sAngles, 360); %get everything within -360:360;
%             tmp(tmp<0) = tmp(tmp<0)+360; %get everything within 0:360
%             tmp(tmp>180) = 360 - tmp(tmp>180); %mirror across x-axis
%             
%             for jj = 1:nBins
%                 tmpCount = sum(tmp>=binLims(jj) & tmp<=binLims(jj+1));
%                 if tmpCount ~= 2
%                     keyboard
%                 end
%             end
            
            s_angles_shuffled = [s_angles_shuffled sAngles]; %#ok<AGROW>
            lam_all_shuff      = [lam_all_shuff lam]; %#ok<AGROW>
        end
    end
    
    %get all angles between 0 and 180
    s_angles_shuffled = rem(s_angles_shuffled, 360); %get everything within -360:360;
    s_angles_shuffled(s_angles_shuffled<0) = s_angles_shuffled(s_angles_shuffled<0)+360; %get everything within 0:360
    s_angles_shuffled(s_angles_shuffled>180) = 360 - s_angles_shuffled(s_angles_shuffled>180); %mirror across x-axis
    
    if params.range90
        s_angles_shuffled(s_angles_shuffled>90) = 180 - s_angles_shuffled(s_angles_shuffled>90); %mirror across y-axis
    end
        
    for ii = 1:nBins
        if ii == 1
            binMeansAllShuffled(kk,ii) = mean(lam_all_shuff(s_angles_shuffled>=binLims(ii) & s_angles_shuffled<=binLims(ii+1)));
            binCountAllShuffled(kk,ii) = sum(s_angles_shuffled>=binLims(ii) & s_angles_shuffled<=binLims(ii+1));
        else
            binMeansAllShuffled(kk,ii) = mean(lam_all_shuff(s_angles_shuffled>binLims(ii) & s_angles_shuffled<=binLims(ii+1)));
            binCountAllShuffled(kk,ii) = sum(s_angles_shuffled>binLims(ii) & s_angles_shuffled<=binLims(ii+1));
        end
        %binMeansPosShuffled(kk,ii) = mean(lam_pos_shuff(s_angles_shuffled>binLims(ii) & s_angles_shuffled<=binLims(ii+1)));
        %binMeansNegShuffled(kk,ii) = mean(lam_neg_shuff(s_angles_shuffled>binLims(ii) & s_angles_shuffled<=binLims(ii+1)));
    end
    
    
    %%%%%%% single-value measures of "peakiness"
    
    % sum of squared residuals
    residAllShuff(kk) = sum((binMeansAllShuffled(kk,:) - lamMean).^2);
        
    % peak bin value
    peakAllShuff(kk) = max(abs(binMeansAllShuffled(kk,:)));
    
    % entropy, in bits
    % separate into positive lambda bins and negative lambda bins
    pAll = [binMeansAllShuffled(kk,:) -binMeansAllShuffled(kk,:)];
    pAll(pAll<0) = 0;
    pAll = pAll/sum(pAll);  %normalize to make them look like probabilities
        
    logPAll = log2(pAll);

    logPAll(pAll==0) = 0; %since limit as p-->0 of plog(p) = 0, replace these log values with zero

    entAllShuff(kk) = -sum(pAll.*logPAll);
end
toc


%% plot results of shuffle analysis

%calculate cumulative means and standard deviations to look for convergence

cumMeansAll = zeros(size(binMeansAllShuffled));
cumStdAll = zeros(size(binMeansAllShuffled));
cumBinCount = zeros(size(binMeansAllShuffled));

cumMeanResid = zeros(1,params.nShuffleReps);
cumStdResid = zeros(1,params.nShuffleReps);


for ii = 1:params.nShuffleReps
    cumMeansAll(ii,:) = mean(binMeansAllShuffled(1:ii,:), 1);
    cumStdAll(ii,:) = std(binMeansAllShuffled(1:ii,:),0, 1);
    cumBinCount(ii,:) = mean(binCountAllShuffled(1:ii,:), 1);

    cumMeanResid(ii) = mean(residAllShuff(1:ii));
    cumStdResid(ii) = std(residAllShuff(1:ii));
        
end

plotCols = lines(size(binMeansAllShuffled,2));
figure
hold on
for ii = 1:size(cumMeansAll,2)
    %plot(cumMeansAll(:,ii), '-', 'color', plotCols(ii,:))
    plot(cumStdAll(2500:end,ii), '--', 'color', plotCols(ii,:))
    %plot(cumBinCount(2500:end,ii)-18, '-', 'color', plotCols(ii,:))
end
xlabel('iteration')
ylabel('cumulative mean/standard deviation for each bin')
title('positive and negative lambdas')


figure
hold on
plot(cumMeanResid, 'k-')
plot(cumStdResid, 'k--')
xlabel('iteration')
ylabel('cumulative mean/standard deviation of sum of squared residuals')
title('all lambdas')



%% plot mean/standard deviations determined by shuffle analysis

%all lambdas
meanCoordsPos = zeros(2,0);
sd1LowPos = zeros(2,0);
sd2LowPos = zeros(2,0);
sd1HighPos = zeros(2,0);
sd2HighPos = zeros(2,0);
meanCoordsNeg = zeros(2,0);
sd1LowNeg = zeros(2,0);
sd2LowNeg = zeros(2,0);
sd1HighNeg = zeros(2,0);
sd2HighNeg = zeros(2,0);

for ii = 1:nBins
    thisMean = cumMeansAll(end,ii);
    thisSd = cumStdAll(end,ii);
    
    arc_th = linspace(binLims(ii), binLims(ii+1), 20);
    
    arc = [cosd(arc_th); sind(arc_th)];
    
    meanCoordsPos = [meanCoordsPos max([0, thisMean])*arc]; 
    sd1LowPos =     [sd1LowPos     max([0, (thisMean-thisSd)])*arc]; 
    sd2LowPos =     [sd2LowPos     max([0, (thisMean-2*thisSd)])*arc]; 
    sd1HighPos =    [sd1HighPos    max([0, (thisMean+thisSd)])*arc]; 
    sd2HighPos =    [sd2HighPos    max([0, (thisMean+2*thisSd)])*arc]; 
    
    meanCoordsNeg = [meanCoordsNeg max([0, -thisMean])*arc]; 
    sd1LowNeg =     [sd1LowNeg     max([0, -(thisMean-thisSd)])*arc]; 
    sd2LowNeg =     [sd2LowNeg     max([0, -(thisMean-2*thisSd)])*arc]; 
    sd1HighNeg =    [sd1HighNeg    max([0, -(thisMean+thisSd)])*arc]; 
    sd2HighNeg =    [sd2HighNeg    max([0, -(thisMean+2*thisSd)])*arc]; 
end

meanAllData = lambdaVAxonPlotterBase(s_angles_all, lam_all, 'binSize', params.binSize,...
    'reflectData', params.reflectData, 'lamLimPlot', params.lamLimPlot);
title('positive and negative lambdas')
hold on

% plot(meanCoordsPos(1,:), meanCoordsPos(2,:), 'k--', 'linewidth', 2)
% plot(sd1LowPos(1,:),  sd1LowPos(2,:),  'k--')
% plot(sd2LowPos(1,:),  sd2LowPos(2,:),  'k:')
% plot(sd1HighPos(1,:), sd1HighPos(2,:), 'k--')
% plot(sd2HighPos(1,:), sd2HighPos(2,:), 'k:')

% plot(meanCoordsNeg(1,:), meanCoordsNeg(2,:), 'r--', 'linewidth', 2)
% plot(sd1LowNeg(1,:),  sd1LowNeg(2,:),  'r--')
% plot(sd2LowNeg(1,:),  sd2LowNeg(2,:),  'r:')
% plot(sd1HighNeg(1,:), sd1HighNeg(2,:), 'r--')
% plot(sd2HighNeg(1,:), sd2HighNeg(2,:), 'r:')


%% plot single-value measures of "peakiness"

%sum of squared residuals measure
residAllData = sum((meanAllData - lamMean).^2);

%peak value measure
peakAllData = max(abs(meanAllData));

%entropy measure (bits)
normAllData = [meanAllData -meanAllData]; % separate into positive lambda bins and negative lambda bins
normAllData(normAllData<0) = 0;
normAllData = normAllData/sum(normAllData);  %normalize to make them look like probabilities

logAllData = log2(normAllData);
logAllData(normAllData==0) = 0;
entAllData = -sum(normAllData.*logAllData);


%calculate means and sigmas for shuffled measures
peakAllShuffMean = mean(peakAllShuff);
peakAllShuffSig  =  std(peakAllShuff);

entAllShuffMean = mean(entAllShuff);
entAllShuffSig  =  std(entAllShuff);

residAllShuffMean = mean(residAllShuff);
residAllShuffSig  =  std(residAllShuff);



figure
subplot(1,3,1)
hold on
xPos = 0.5;
%all lambdas, sum of squared residuals measure
plot(xPos*[1 1], residAllShuffMean + residAllShuffSig*[-1 1], 'r-', 'linewidth', 3)
plot(xPos*[1 1], residAllShuffMean + 2*residAllShuffSig*[-1 1], 'r-', 'linewidth', 1)
plot(xPos, residAllShuffMean, 'ro', 'MarkerFaceColor', [1 1 1])
plot(xPos, residAllData, 'k.', 'MarkerFaceColor', [0 0 0])

subplot(1,3,2)
hold on

%all lambdas, peak measure
xPos = 0;
plot(xPos*[1 1], peakAllShuffMean + peakAllShuffSig*[-1 1], 'r-', 'linewidth', 3)
plot(xPos*[1 1], peakAllShuffMean + 2*peakAllShuffSig*[-1 1], 'r-', 'linewidth', 1)
plot(xPos, peakAllShuffMean, 'ro', 'MarkerFaceColor', [1 1 1])
plot(xPos, peakAllData, 'k.', 'MarkerFaceColor', [0 0 0])

xlabel('all lambdas')
set(gca, 'xLim', [-1 1])
title('peak value measure')


subplot(1,3,3)
hold on

%positive and negative lambdas, entropy
xPos = 0;
plot(xPos*[1 1], entAllShuffMean + entAllShuffSig*[-1 1], 'r-', 'linewidth', 3)
plot(xPos*[1 1], entAllShuffMean + 2*entAllShuffSig*[-1 1], 'r-', 'linewidth', 1)
plot(xPos, entAllShuffMean, 'ro', 'MarkerFaceColor', [1 1 1])
plot(xPos, entAllData, 'k.', 'MarkerFaceColor', [0 0 0])

xlabel('all lambdas')
set(gca, 'xLim', [-1 1])
title('entropy measure')




%% plot distributions from random rotation analysis


xMax = (1/100)*ceil(max([residAllShuff; residAllData])*100);
binEdges = 0:0.0002:xMax;
binCounts = histc(residAllShuff, binEdges);

%plot histogram
patchVals = [binEdges(1); 0];

for ii = 1:length(binEdges)-1
    patchVals = [patchVals, [binEdges(ii) binEdges(ii+1); binCounts(ii) binCounts(ii)]];
end
patchVals = [patchVals, [binEdges(end); 0]];

figure('position', [200 200 400 200]); hold on
patch(patchVals(1,:), patchVals(2,:), 'r')
plot(residAllData, 100, 'k*')

nAboveObs = sum(residAllShuff>residAllData);

xlabel('anisotropy index')
title(['random >  observed: ' num2str(nAboveObs) ' of ' num2str(params.nShuffleReps)])
set(gca, 'xlim', [0 xMax])



