function [re_master_datarun, re_slave_datarun] = map(master_datarun, slave_datarun, varargin)
% map_gdf       wrapper function for EI based mapping
%                   note: work in progress, must be extended
%                   just takes best map from map_ei and organises cell_types 
%
% usage:   [datarun_1, datarun_2] = map_gdf(datarun_1, datarun_2, varargin)
%
% arguments:    master_datarun - datarun structure for the master data run
%                slave_datarun - datarun structure for the slave data run 
%
% outputs:  re_master_datarun - datarun structure with increased datarun     
%
% optional fields with defaults
%   
%   master_cell_type        'all'     cell type identifier for master datarun
%                                     (see get_cell_indices for more information)
%   slave_cell_type         'all'     cell type identifier for slave datarun
%   electrode_threshold       5       electrodes with signals below electrode_threshold are set to 0
%   significant_electrodes    10      there must be at least this number of electrodes with a significant signal
%   corr_threshold           0.85     only cell pairs with correlation coefficients above this 
%                                     threshold are reported
%   space_only              true      use only the spatial information in the EI
%   verbose                 false     logical for producing verbose output
%   troubleshoot             false     logical for printing a cell by cell
%                                     account of the mapping to the command
%                                     window
%
% 2008-10
%   GDF modification of function made by MG


% SET UP OPTIONAL ARGUMENTS

p = inputParser;
p.addRequired('master_datarun', @isstruct)
p.addRequired('slave_datarun', @isstruct)
p.addParamValue('master_cell_type', 'all')
p.addParamValue('slave_cell_type', 'all')
p.addParamValue('electrode_threshold', 5, @isnumeric)
p.addParamValue('significant_electrodes', 10, @isnumeric)
p.addParamValue('corr_threshold', 0.85, @isnumeric)
p.addParamValue('space_only', true, @islogical);
p.addParamValue('verbose', false, @islogical)
p.addParamValue('troubleshoot', false, @islogical)
p.parse(master_datarun, slave_datarun, varargin{:});

verbose = p.Results.verbose;
master_cell_type = p.Results.master_cell_type;
slave_cell_type = p.Results.slave_cell_type;
troubleshoot = p.Results.troubleshoot;

% generate a struct to pass varagin from wrapper to map_ei
map_ei_params = rmfield(p.Results, {'master_datarun', 'slave_datarun'});

cell_list_map = map_ei(master_datarun, slave_datarun, map_ei_params);

re_master_datarun = master_datarun;
re_slave_datarun = slave_datarun;

for i = 1:length(master_cell_type)
    
    if ~isnumeric(master_cell_type{i})
        for ii = 1:length(master_datarun.cell_types)
            temp = isequal(master_cell_type{i}, master_datarun.cell_types{ii}.name);
        end
        cell_types_nr = find(temp);
    else
        cell_types_nr = master_cell_type{i};
    end

    cell_types_1.name = master_datarun.cell_types{cell_types_nr}.name;
    cell_types_2.name = master_datarun.cell_types{cell_types_nr}.name;
    cell_types_1.cell_ids = [];
    cell_types_2.cell_ids = [];    
        
    for ii = 1:length(master_datarun.cell_types{cell_types_nr}.cell_ids)
        t = find(master_datarun.cell_ids == master_datarun.cell_types{cell_types_nr}.cell_ids(ii));
        if ~isempty(cell_list_map{t})            
            cell_types_1.cell_ids = [cell_types_1.cell_ids master_datarun.cell_types{cell_types_nr}.cell_ids(ii)]; 
            cell_types_2.cell_ids = [cell_types_2.cell_ids cell_list_map{t}(1)]; 
        end
    end

        
    if verbose
        fprintf('mapping: %s \n', cell_types_1.name)
    	fprintf('Of %3.0f cells, %3.0f pairs were found in the slave \n',length(master_datarun.cell_types{cell_types_nr}.cell_ids),length(cell_types_1.cell_ids));
    end
  
    
    re_master_datarun.cell_types{cell_types_nr} = cell_types_1;
    re_slave_datarun.cell_types{cell_types_nr} = cell_types_2; 
    re_slave_datarun.duplicates = cell_list_map;
    
end

