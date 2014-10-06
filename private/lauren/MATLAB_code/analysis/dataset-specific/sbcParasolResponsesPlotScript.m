%% script to plot response curves of midgets and parasols from 2008-08-27-2-data002

clear all
close all

%analysis parameters
recalcAll = false; %redo all response curve fits and standard deviation calculations
nBootstrapReps = 100;

maxSuccRateCutoff = 0.4;

if 0
    analysisPath = '/snle/lab/Experiments/Array/Analysis/2011-10-25-4/data002/';
    fileNameSuffix = '_w50'; %use when _w50 or _w100 is included in filenames
    
    pathToEi = '/snle/lab/Experiments/Array/Analysis/2011-10-25-4/data000/data000.ei';
    %datarunPath = '2011-10-25-4/data000/data000';
    
    cellInfo = cell_list_2011_10_25_4();
    %%% criteria for included cells: max (negative) signal must either
    % (1) not fall on an edge electrode OR
    % (2) be >= cutoffMult of the mean max signal of cells within class that meet (1)

    [cellInfo cellInfoOffArray] = removeSmallSigEdgeCells(cellInfo, pathToEi, 'cutoffMult', 0.5);

    
else
    analysisPath = '/snle/lab/Experiments/Array/Analysis/2012-01-27-3/data002/';
    fileNameSuffix = '_w100'; %use when _w50 or _w100 is included in filenames
    
    %pathToEi = '/snle/lab/Experiments/Array/Analysis/2012-01-27-3/data000/data000.ei';
    %datarunPath = '2012-01-27-3/data000/data000';
    
    cellInfo = cell_list_2012_01_27_3();
    
    %%% criteria for included cells: max (negative) signal must either
    % (1) not fall on an edge electrode OR
    % (2) be >= cutoffMult of the mean max signal of cells within class that meet (1)
    
    %do parasols and sbcs separately because they use different ei files
    iOP = 0; iSBC = 0;
    for ii = 1:length(cellInfo)
        if strcmpi(cellInfo(ii).type, 'onPar')
            iOP = iOP+1;
            cellInfoPar(iOP) = cellInfo(ii);
        elseif strcmpi(cellInfo(ii).type, 'sbc')
            iSBC = iSBC+1;
            cellInfoSBC(iSBC) = cellInfo(ii);
        else
            error('unexpected cell type')
        end
    end
    [cellInfoPar cellInfoParOffArray] = removeSmallSigEdgeCells(cellInfoPar, cellInfoPar(1).eiPath, 'cutoffMult', 0.5);
    [cellInfoSBC cellInfoSBCOffArray] = removeSmallSigEdgeCells(cellInfoSBC, cellInfoSBC(1).eiPath, 'cutoffMult', 0.5);

    %recombine structs
    cellInfo = [cellInfoSBC cellInfoPar];
    cellInfoOffArray = [cellInfoSBCOffArray cellInfoParOffArray];
end

arraySpacing = 60; %pitch of electrodes in microns
cutoffMult = 0.5; %determines which cells will be included in analysis

%colors!
blue = [50 70 247]/255;
rust = [.8 .05 0.05];
grass = [90 156 0]/255;
salmon = [255 124 59]/255;

cellTypes = {'onPar', 'sbc'};

colors(1,:) = 1-(1-blue)*0.7; %ON parasol
colors(2,:) = grass;


%% separate out unstimulated cells into a list (remove from cellInfo struct)

count = 0;
for ii = length(cellInfo):-1:1
    if isempty(cellInfo(ii).stimElec)
        count = count + 1;
        cellInfoNoStim(count) = cellInfo(ii);
        cellInfo(ii) = [];
    end
end


%% getting data from elecResp files

nCells = length(cellInfo);

for ii = 1:nCells
    load([analysisPath filesep 'elecResp_n' num2str(cellInfo(ii).id) '_p' num2str(cellInfo(ii).stimElec) fileNameSuffix '.mat'])
    elecResp = checkForUnfinishedAnalysis(elecResp, nBootstrapReps, 'recalcAll', recalcAll);
    save([analysisPath filesep 'elecResp_n' num2str(cellInfo(ii).id) '_p' num2str(cellInfo(ii).stimElec) fileNameSuffix '.mat'], 'elecResp')
    
    %get response rates
    cellInfo(ii).data = zeros(3, length(elecResp.stimInfo.stimAmps));
    cellInfo(ii).data(1,:) = abs(elecResp.stimInfo.stimAmps);
    cellInfo(ii).data(2,:) = elecResp.analysis.successRates;
    cellInfo(ii).data(3,:) = elecResp.stimInfo.nPulses;
    
    %remove unanalyzed data
    for jj = length(elecResp.stimInfo.stimAmps):-1:1
        if isempty(elecResp.analysis.type{jj})
            cellInfo(ii).data(:,jj) = [];
        end
    end
    
    %get fit parameters
    cellInfo(ii).params = elecResp.analysis.erfParams;
    cellInfo(ii).SD = 1/(sqrt(2)*cellInfo(ii).params(1)); %standard deviation of cumulative Gaussian fit, in µA
    cellInfo(ii).thresh = elecResp.analysis.threshold; % should equal -cellInfo(ii).params(2)/cellInfo(ii).params(1)
end

%% apply maxSuccRateCutOff

count = 0;
for ii = length(cellInfo):-1:1
    if max(cellInfo(ii).data(2,:)) < maxSuccRateCutoff

        disp(['cell ' num2str(cellInfo(ii).id) ' didn''t meet success rate requirement of at least ' num2str(maxSuccRateCutoff)])
        
        %add cell to "cellInfoNoStim"
        ind = length(cellInfoNoStim)+1;
        cellInfoNoStim(ind).id = cellInfo(ii).id;
        cellInfoNoStim(ind).stimElec = cellInfo(ii).stimElec;
        cellInfoNoStim(ind).verMin = cellInfo(ii).verMin;
        cellInfoNoStim(ind).type = cellInfo(ii).type;
        
        %remove cell from "cellInfoNoStim"
        cellInfo(ii) = [];
    end
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

%% summary plot 


params.threshBarLim = [0 250];
params.PW = 100;
params.chargeLim = 1.9075*params.PW;
params.safetyLim = 135.0;

params.yBarPos = [0.55 0.2];
params.OS = 0.12;


figure('position', [100 100 800 220], 'color', [1 1 1])
axes('position', [0.1 0.1 0.8 0.8])
make_1D_thresh_plot(cellInfo, colors, params, cellTypes)


