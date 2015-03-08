clear
datapath='/Volumes/Analysis/nora/NSEM/GLM_Output/fixedSP_rk1_linear_MU_PS_CP_p8IDp8/standardparams/';
exp_names=['2012-08-09-3/';'2012-09-27-3/';'2013-08-19-6/';'2013-10-10-0/'];
fittypepath{2}='NSEM_mapPRJ/';
fittypepath{1}='WN_mapPRJ/';

data = zeros(10000, 4);
data_start = 1;

for exp=1
    for fittype=1:2
        
        % Get file list
        matfiles=dir([datapath fittypepath{fittype} exp_names(exp,:) '*.mat']);
        n_files=length(matfiles);
        
        % Initialize matrix
        cell_number{fittype} = zeros(n_files,1);
        
        % get all cell ids for a given experiment
        for file=1:n_files
            load([datapath fittypepath{fittype} exp_names(exp,:) matfiles(file).name]);
            cell_number{fittype}(file) = fittedGLM.cellinfo.cid;
        end
        
    end
    
    % only include cells in both WN and NSEM runs
    cells = intersect(cell_number{1}, cell_number{2}, 'stable');
    n_cells = length(cells);

    % initialize
    coupling = zeros(2,n_cells, n_cells);
    types=zeros(n_cells,1);
    cell_loc = zeros(n_cells, 2);
    
    for fittype=1:2
        
        % Get file list
        matfiles=dir([datapath fittypepath{fittype} exp_names(exp,:) '*.mat']);
        n_files=length(matfiles);
        
        % Collect coupling info from files
        for file=1:n_files
            load([datapath fittypepath{fittype} exp_names(exp,:) matfiles(file).name]);
            idx = find(cells == fittedGLM.cellinfo.cid);
            if ~isempty(idx)
                
                for pair=1:12
                    CP_idx=find(cells == fittedGLM.cellinfo.pairs(pair));
                    coupling(fittype, idx, CP_idx) = max(fittedGLM.linearfilters.Coupling.Filter{pair});
                end
                
                if strcmp(matfiles(file).name(2), 'N')
                    types(idx)=1/2;
                else
                    types(idx)=-1/2;
                end
                
                cell_loc(idx, 1) = fittedGLM.cellinfo.slave_centercoord.x_coord;
                cell_loc(idx, 2) = fittedGLM.cellinfo.slave_centercoord.y_coord;
                
            end
        end
        
    end
    
    distances = squareform(pdist(cell_loc));
    CP_WN=coupling(1,:,:);
    CP_NSEM=coupling(2,:,:);
    
    data_idx = find(CP_WN ~= 0);
    test = find (CP_NSEM ~= 0);
    if sum(test == dataidx) ~= 0
        warn('shit')
    end
    Data(data_start:(data_start+length(data_idx)),1) = distances(data_idx);
    Data(data_start:(data_start+length(data_idx)),2) = type(data_idx);
    Data(data_start:(data_start+length(data_idx)),3) = CP_WN(data_idx);
    Data(data_start:(data_start+length(data_idx)),4) = CP_NSEM(data_idx);
    data_start = data_start+length(data_idx)+1;

end
