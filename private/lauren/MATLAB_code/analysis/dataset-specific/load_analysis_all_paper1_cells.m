function cellInfoByType = load_analysis_all_paper1_cells(varargin)

% load basic analysis information for all cells that go into
% midget/SBC/parasol manuscript
%
% useful for plotting basic stuff about all cells such as latencies and
% stim elec positions
%
% ***note that this does not load up cells from cell_list_low_freq_stim
% (cells used to check for low-latency responses)
%


p = inputParser;

p.addParamValue('recalcAll', false, @islogical)
p.addParamValue('nBootstrapReps', 100, @isnumeric)

%minimum response probability in an analyzed movie for a cell to be
%included as "stimulated"
p.addParamValue('maxSuccRateCutoff', 0.4, @isnumeric)

%whether to remove cells that are determined to be positioned off the array from datasets
%involving all cells of a given type in a given prep
p.addParamValue('removeOffArrayCells', true, @islogical)

%if removeOffArrayCells == true, use this value as the cutoff ei signal
%for edge cells (see removeSmallSigEdgeCells)
p.addParamValue('cutoffMult', 0.5, @isnumeric)

%whether to exclude sbc data that was recorded on 30 micron arrays
p.addParamValue('exclude30ArrayData', true, @islogical)

%which of the datasets to include (see details below)
p.addParamValue('includeDataset', true(1,7), @islogical)

p.addParamValue('cellTypes', {'onPar', 'offPar', 'onMidg', 'offMidg', 'sbc'}, @iscell)

p.parse(varargin{:})

params = p.Results;


%% dataset details

analysisPath{1} = '/snle/lab/Experiments/Array/Analysis/2008-08-27-2/data002/';
fileNameSuffix{1} = ''; %use when _w50 or _w100 is included in filenames
PW{1} = 100;
pathToEi{1} = '/snle/lab/Experiments/Array/Analysis/2008-08-27-2/data001-lh/data001-lh.ei';
cellInfo{1} = cell_list_2008_08_27_2();

analysisPath{2} = '/snle/lab/Experiments/Array/Analysis/2011-01-11-0/data029/';
fileNameSuffix{2} = '_w50'; %use when _w50 or _w100 is included in filenames
PW{2} = 50;
pathToEi{2} = '/snle/lab/Experiments/Array/Analysis/2011-01-11-0/data030-lh/data030-lh.ei';
cellInfo{2} = cell_list_2011_01_11_0();

analysisPath{3} = '/snle/lab/Experiments/Array/Analysis/2011-10-25-4/data002/';
fileNameSuffix{3} = '_w50'; %use when _w50 or _w100 is included in filenames
PW{3} = 50;
pathToEi{3} = '/snle/lab/Experiments/Array/Analysis/2011-10-25-4/data000/data000.ei';
cellInfo{3} = cell_list_2011_10_25_4();

analysisPath{4} = '/snle/lab/Experiments/Array/Analysis/2012-01-27-3/data002/';
PW{4} = 100;
fileNameSuffix{4} = '_w100';
pathToEi{4} = []; %different for different cell types
cellInfo{4} = cell_list_2012_01_27_3();

%collection of SBCs from various preps
analysisPath{5} = '/snle/lab/Experiments/Array/Analysis/';
prefPW = 100;
cellInfo{5} = cell_list_sbc(params.exclude30ArrayData);

%cells plotted in overview figure
analysisPath{6} = '/snle/lab/Experiments/Array/Analysis/';
cellInfo{6} = cell_list_overview_fig();

%cells investigated for long-latency responses
[cellInfoTmp1 cellInfoTmp2] = cell_list_low_freq_stim();
cellInfoTmp1 = rmfield(cellInfoTmp1, {'path', 'movieTimePath', 'elecRespPath', 'WN'});
cellInfoTmp2 = rmfield(cellInfoTmp2, {'splitElecRespPath', 'originalElecRespPath', 'repStartTimesPath', 'movieTimePath'});
cellInfo{7} = [cellInfoTmp1 cellInfoTmp2];
clear cellInfoTmp1 cellInfoTmp2

%% apply exclusions

%%% separate out off-array cells for datasets corresponding to single
%%% preparations
if params.removeOffArrayCells
    for iDS = 1:3
        [cellInfo{iDS} cellInfoOffArray{iDS}] = removeSmallSigEdgeCells(cellInfo{iDS}, pathToEi{iDS}, 'cutOffMult', params.cutoffMult);
    end
    
    %special case: 4th dataset (different ei files for different cell types)
    iDS = 4;
    iOP = 0; iSBC = 0;
    for ii = 1:length(cellInfo{iDS})
        if strcmpi(cellInfo{iDS}(ii).type, 'onPar')
            iOP = iOP+1;
            cellInfoPar(iOP) = cellInfo{iDS}(ii);
        elseif strcmpi(cellInfo{iDS}(ii).type, 'sbc')
            iSBC = iSBC+1;
            cellInfoSBC(iSBC) = cellInfo{iDS}(ii);
        else
            error('unexpected cell type')
        end
    end
    [cellInfoPar cellInfoParOffArray] = removeSmallSigEdgeCells(cellInfoPar, cellInfoPar(1).eiPath, 'cutoffMult', 0.5);
    [cellInfoSBC cellInfoSBCOffArray] = removeSmallSigEdgeCells(cellInfoSBC, cellInfoSBC(1).eiPath, 'cutoffMult', 0.5);

    %recombine structs
    cellInfo{iDS} = [cellInfoSBC cellInfoPar];
    cellInfoOffArray{iDS} = [cellInfoSBCOffArray cellInfoParOffArray];
end



for iDS = 1:7
    %make sure all cells have been analyzed
    for ii = 1:length(cellInfo{iDS})
        if cellInfo{iDS}(ii).stimElec == 0 %indicates unfinished analysis
            error(['unfinished analysis for cell ' num2str(cellInfo{iDS}(ii).id) ' from ' cellInfo{iDS}(ii).analysisPath])
        end
    end
    
    %%% remove unstimulated cells
    for ii = length(cellInfo{iDS}):-1:1
        if isempty(cellInfo{iDS}(ii).stimElec) %indicates that no stimulation by any electrode was detected
            cellInfo{iDS}(ii) = [];
        end
    end
end


%% get analysis

for iDS = 1:7
    if params.includeDataset(iDS)
        % load analysis from elecResp files
        for ii = 1:length(cellInfo{iDS})
            
            if any([1 2 3 4] == iDS)
                fullPath = [analysisPath{iDS} filesep 'elecResp_n' num2str(cellInfo{iDS}(ii).id)...
                    '_p' num2str(cellInfo{iDS}(ii).stimElec) fileNameSuffix{iDS} '.mat'];
                cellInfo{iDS}(ii).PW = PW{iDS};
            elseif iDS == 5 %mixture of SBCs from different preps
                if length(cellInfo{iDS}(ii).PW) > 1
                    %%% use other PW if preferred PW doesn't fulfill maxSuccRateCutoff %%%
                    tmp = load([analysisPath{iDS} filesep cellInfo{iDS}(ii).analysisPath filesep...
                        'elecResp_n' num2str(cellInfo{iDS}(ii).id) '_p' num2str(cellInfo{iDS}(ii).stimElec) '_w' num2str(prefPW) '.mat']);
                    
                    %determine maximum measured response rate
                    succRatesTmp = tmp.elecResp.analysis.successRates;
                    for jj = length(tmp.elecResp.stimInfo.stimAmps):-1:1
                        if isempty(tmp.elecResp.analysis.type{jj})
                            succRatesTmp(jj) = [];
                        end
                    end
                    
                    if max(succRatesTmp) >= params.maxSuccRateCutoff
                        cellInfo{iDS}(ii).PW = prefPW;
                    else
                        cellInfo{iDS}(ii).PW = cellInfo{iDS}(ii).PW(cellInfo{iDS}(ii).PW~=prefPW);
                        disp(['Preferred pulse width doesn''t satisfy maxSuccRateCutoff criterion for cell ' num2str(cellInfo{iDS}(ii).id) ': using other pulse width (' num2str(cellInfo{iDS}(ii).PW) ')'])
                    end
                end
                
                %load elecResp (look for elecResp with specified PW value first)
                if exist([analysisPath{iDS} filesep cellInfo{iDS}(ii).analysisPath filesep...
                        'elecResp_n' num2str(cellInfo{iDS}(ii).id) '_p' num2str(cellInfo{iDS}(ii).stimElec) '_w' num2str(cellInfo{iDS}(ii).PW) '.mat'], 'file')
                    fullPath = [analysisPath{iDS} filesep cellInfo{iDS}(ii).analysisPath filesep...
                        'elecResp_n' num2str(cellInfo{iDS}(ii).id) '_p' num2str(cellInfo{iDS}(ii).stimElec) '_w' num2str(cellInfo{iDS}(ii).PW) '.mat'];
                elseif exist([analysisPath{iDS} filesep cellInfo{iDS}(ii).analysisPath filesep...
                        'elecResp_n' num2str(cellInfo{iDS}(ii).id) '_p' num2str(cellInfo{iDS}(ii).stimElec) '.mat'], 'file');
                    fullPath = [analysisPath{iDS} filesep cellInfo{iDS}(ii).analysisPath filesep...
                        'elecResp_n' num2str(cellInfo{iDS}(ii).id) '_p' num2str(cellInfo{iDS}(ii).stimElec) '.mat'];
                    disp(['no PW specifier in file elecResp_n' num2str(cellInfo{iDS}(ii).id) '_p' num2str(cellInfo{iDS}(ii).stimElec)])
                else
                    error(['could not find file elecResp_n' num2str(cellInfo{iDS}(ii).id) '_p' num2str(cellInfo{iDS}(ii).stimElec)])
                end
            elseif iDS == 6
                fullPath = [cellInfo{iDS}(ii).pathToAnalysis filesep 'elecResp_n' num2str(cellInfo{iDS}(ii).id)...
                    '_p' num2str(cellInfo{iDS}(ii).stimElec) cellInfo{iDS}(ii).suffix '.mat'];
            elseif iDS == 7
                fullPath = cellInfo{iDS}(ii).scanElecRespPath;
            else %iDS == 8
                
            end
            
            load(fullPath);
            elecResp = checkForUnfinishedAnalysis(elecResp, params.nBootstrapReps, 'recalcAll', params.recalcAll, 'checkForAnalysisGaps', false);
            save(fullPath, 'elecResp')
            
            cellInfo{iDS}(ii).elecRespPath = fullPath;
            clear fullPath
            
            cellInfo{iDS}(ii).thresh = elecResp.analysis.threshold; % should equal -cellInfo{iDS}(ii).params(2)/cellInfo{iDS}(ii).params(1)
            cellInfo{iDS}(ii).maxSuccRate = max(elecResp.analysis.successRates); %#ok<*AGROW>
            cellInfo{iDS}(ii).ei = elecResp.cells.mainEI;
            cellInfo{iDS}(ii).shortName = elecResp.names.rrs_short_name;
            cellInfo{iDS}(ii).experiment = elecResp.names.experiment;
            cellInfo{iDS}(ii).elecRespID = elecResp.cells.main;
            
            tmp = strfind(elecResp.names.rrs_ei_path, filesep);
            cellInfo{iDS}(ii).eiFile = elecResp.names.rrs_ei_path(tmp(end)+1:end);
            
            clear elecResp
            
        end
    else
        cellInfo{iDS} = [];
    end
end

%% remove duplicate cells that show up in more than one dataset (keeps cell
% info from first dataset containing cell)

for ii = 2:length(cellInfo)
    dupCells = false(1, length(cellInfo{ii}));
    for jj = 1:length(cellInfo{ii})
        for kk = 1:ii-1 %just check preceding datasets
            if ~isempty(cellInfo{kk})
                %matches = cellfun(@(x)strcmp(x,cellInfo{ii}(jj).shortName), {cellInfo{kk}.shortName});
                %check for match of piece and cell (doesn't require matched stim dataset, PW, or stim electrode)
                expMatches = cellfun(@(x)strcmp(x,cellInfo{ii}(jj).experiment), {cellInfo{kk}.experiment});
                eiFileMatches = cellfun(@(x)strcmp(x,cellInfo{ii}(jj).eiFile), {cellInfo{kk}.eiFile});
                idMatches = cellInfo{ii}(jj).elecRespID == [cellInfo{kk}.elecRespID];
                typeMatches = cellfun(@(x)strcmpi(x,cellInfo{ii}(jj).type), {cellInfo{kk}.type});
                
                if any(expMatches & eiFileMatches & idMatches) %has to be the same cell
                    dupCells(jj) = true;
                    disp(['removed cell ' cellInfo{ii}(jj).shortName ' from dataset ' num2str(ii) ' because same cell is already in dataset ' num2str(kk)])
                    
                    stimElec1 = cellInfo{ii}(jj).stimElec;
                    stimElec2 = cellInfo{kk}(expMatches & eiFileMatches & idMatches).stimElec;
                    if stimElec1 ~= stimElec2
                        disp(['*** WARNING: same cell appears in 2 datasets with different stim elecs (' num2str(stimElec1) ' and ' num2str(stimElec2) ')'])
                    end
                elseif any(expMatches & typeMatches & ~eiFileMatches) %can't determine whether cell is the same -- need to check manually
                    disp(['may be a duplicate cell in dataset ' num2str(ii) ': ',...
                        'need to check whether neuron ' num2str(cellInfo{ii}(jj).elecRespID) ' in ' cellInfo{ii}(jj).experiment ' ' cellInfo{ii}(jj).eiFile...
                        ' is already represented in dataset ' num2str(kk)])
                end
            end
        end
    end
    cellInfo{ii}(dupCells) = [];
end


%% apply maxSuccRateCutOff

for jj = 1:length(cellInfo)
    for ii = length(cellInfo{jj}):-1:1
        if cellInfo{jj}(ii).maxSuccRate < params.maxSuccRateCutoff
            disp(['cell ' num2str(cellInfo{jj}(ii).id) ' didn''t meet success rate requirement of at least ' num2str(params.maxSuccRateCutoff)])
            cellInfo{jj}(ii) = [];
        end
    end
end

%% get rid of inconsistent fields

fieldsToKeep = {'id', 'stimElec', 'type', 'axonElecs', 'PW', 'elecRespPath', 'thresh', 'maxSuccRate', 'ei', 'shortName', 'experiment', 'elecRespID', 'eiFile'};

for ii = 1:length(cellInfo)
    if ~isempty(cellInfo{ii})
        fNames = fieldnames(cellInfo{ii});
        for jj = 1:length(fNames)
            if ~any(strcmp(fNames{jj}, fieldsToKeep))
                cellInfo{ii} = rmfield(cellInfo{ii}, fNames{jj});
            end
        end
        
        if ~isfield(cellInfo{ii}, 'axonElecs')
            cellInfo{ii}(1).axonElecs = [];
        end
        cellInfo{ii} = orderfields(cellInfo{ii});
    end
end

%% collect info by cell type

cellInfoByType = cell(5,1);
for ii = 1:length(cellInfo)
    for jj = 1:length(cellInfo{ii})
        cellTypeInd = find(cellfun(@(x)strcmpi(x,cellInfo{ii}(jj).type), params.cellTypes));
        
        if length(cellTypeInd)==1
            cellInfoByType{cellTypeInd}(length(cellInfoByType{cellTypeInd})+1) = cellInfo{ii}(jj);
        else
            error('invalid cell type')
        end
    end
end





end

