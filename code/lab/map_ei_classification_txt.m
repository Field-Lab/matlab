function map_ei_classification_txt(master_datarun, slave_datarun, varargin)
% map       wrapper function for EI based mapping
%
% try to map as many cells (master_cell_type) as possible, but exclude duplicates 
%
% can handel only two datasets momentarily
%
%   greschner modification from map.m

%example
if 0
    clear datarun
    if 1
        datarun{1}.names.rrs_params_path='/Data.noindex/Greschner/2010-03-05-2/data026/data026.params';
        datarun{1}.names.rrs_neurons_path='/Data.noindex/Greschner/2010-03-05-2/data026/data026.neurons';
        datarun{1}.names.rrs_ei_path='/Data.noindex/Greschner/2010-03-05-2/data026/data026.ei';
        datarun{1}.piece.array_id=1054;
        datarun{2}.names.rrs_params_path='/Data.noindex/Greschner/2010-03-05-2/data025/data025.params';
        datarun{2}.names.rrs_neurons_path='/Data.noindex/Greschner/2010-03-05-2/data025/data025.neurons';
        datarun{2}.names.rrs_ei_path='/Data.noindex/Greschner/2010-03-05-2/data025/data025.ei';
        datarun{2}.piece.array_id=1054;

        save_path='/Data.noindex/Greschner/2010-03-05-2/data025/';
    end    

    opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_ei',1);
    datarun=load_data(datarun,opt);

    % try to map as many cells (master_cell_type) as possible, but exclude duplicates 
    map_ei_classification_txt(datarun{1}, datarun{2},'master_cell_type',{1 2 3 4 5}, 'classification_path', save_path);
end






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
p.addParamValue('verbose', true, @islogical)
p.addParamValue('troubleshoot', true, @islogical)
p.addParamValue('classification_path', '')
p.addParamValue('classification_file_name', 'classification_ei_map-from-')
p.addParamValue('map_file_name', 'ei_map-from-')
p.parse(master_datarun, slave_datarun, varargin{:});

verbose = p.Results.verbose;
master_cell_type = p.Results.master_cell_type;
troubleshoot = p.Results.troubleshoot;

% generate a struct to pass varagin from wrapper to map_ei
map_ei_params = rmfield(p.Results, {'classification_path', 'classification_file_name', 'map_file_name'});
cell_list_map = map_ei(master_datarun, slave_datarun, map_ei_params);

re_master_datarun = master_datarun;
re_slave_datarun = slave_datarun;

res=cell(length(slave_datarun.cell_ids),1);
master_cell_ids=[];
slave_cell_ids=[];

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
        fprintf('\nmapping: %s \n', cell_types_1.name)
    	fprintf('of %3.0f cells, %3.0f pairs were found in the slave \n',length(master_datarun.cell_types{cell_types_nr}.cell_ids),length(cell_types_1.cell_ids));
    end
  

    for i=1:length(cell_types_2.cell_ids)
        t=find(slave_datarun.cell_ids == cell_types_2.cell_ids(i));
        res{t}=cell_types_2.name;
    end
    
    master_cell_ids=[master_cell_ids cell_types_1.cell_ids];
    slave_cell_ids=[slave_cell_ids cell_types_2.cell_ids];
    
end


%write file
    [fid1,msg]=fopen([p.Results.classification_path p.Results.classification_file_name master_datarun.names.short_name '.txt'], 'w');
    for i=1:length(slave_datarun.cell_ids)
        if ~isempty(res{i})
            fprintf(fid1,'%d  All/%s\n',slave_datarun.cell_ids(i),res{i});
        else
            fprintf(fid1,'%d  All\n',slave_datarun.cell_ids(i));
        end
    end
    fclose(fid1);
    
    [fid2,msg]=fopen([p.Results.classification_path p.Results.map_file_name master_datarun.names.short_name '.txt'], 'w');
    for i=1:length(master_cell_ids)
    	fprintf(fid2,'%d\t%d\n',master_cell_ids(i),slave_cell_ids(i));
    end
    fclose(fid2);

   if verbose
        fprintf('\nwrote %s \n', [p.Results.classification_path p.Results.classification_file_name master_datarun.names.short_name '.txt'])
        fprintf('\nwrote %s \n', [p.Results.classification_path p.Results.map_file_name master_datarun.names.short_name '.txt'])
   end






