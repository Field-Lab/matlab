function [cell_list_map, failed_cells, output_matrix] = map_ei(master_datarun, slave_datarun, varargin)
% 
% map_ei
%
% map cells from master_datarun to slave_datarun by the similarity of their
% EIs.
%
% usage: [cell_list_map, failed_cells] = map_ei(master_datarun,slave_datarun,varargin)
%
% arguments:    master_datarun - datarun structure of the master data run
%                slave_datarun - datarun structure of the slave data run
%
% outputs:       cell_list_map - list of indices associating cells in the 
%                                slave datarun with those in the master datarun 
%                 failed_cells - a list of cell IDs from the master without
%                                an identified pair in the slave
%
% optional fields with defaults
%   
%   master_cell_type        'all'     cell type identifier for master datarun
%                                     (see get_cell_indices for more information)
%   slave_cell_type         'all'     cell type identifier for slave datarun
%   electrode_threshold       5       electrodes with signals below electrode_threshold are set to 0
%   significant_electrodes    10      there must be at least this number of electrodes with a significant signal
%   corr_threshold           0.95     only cell pairs with correlation coefficients above this 
%                                     threshold are reported
%   space_only              true      use only the spatial information in the EI
%   verbose                 false     logical for producing verbose output
%   troubleshoot            false     logical for printing a cell by cell
%                                     account of the mapping to the command window
%
% 2008-10 GDF cleanup of function written by MG
% 2010-03 greschner: clean up,  if largest_corr==largest_corr_col

% SET UP OPTIONAL ARGUMENTS
p = inputParser;
p.addRequired('master_datarun', @isstruct);
p.addRequired('slave_datarun',  @isstruct);
p.addParamValue('significant_electrodes',   10,     @isnumeric);
p.addParamValue('master_cell_type',         'all'             );
p.addParamValue('slave_cell_type',          'all'             );
p.addParamValue('electrode_threshold',      5,      @isnumeric);
p.addParamValue('corr_threshold',           0.95,   @isnumeric);
p.addParamValue('verbose',                  false,  @islogical);
p.addParamValue('space_only',               true,   @islogical);
p.addParamValue('troubleshoot',             false,  @islogical);
p.parse(master_datarun, slave_datarun, varargin{:});

verbose = p.Results.verbose;
master_cell_type = p.Results.master_cell_type;
slave_cell_type = p.Results.slave_cell_type;
electrode_threshold = p.Results.electrode_threshold;
corr_threshold = p.Results.corr_threshold;
space_only = p.Results.space_only;
significant_electrodes = p.Results.significant_electrodes;
troubleshoot = p.Results.troubleshoot;
        
% get cell numbers
master_indices = get_cell_indices(master_datarun, master_cell_type);
slave_indices = get_cell_indices(slave_datarun, slave_cell_type);

if ~space_only
    len = min(size(master_datarun.eis{master_indices(1)},2),size(slave_datarun.eis{slave_indices(1)},2));

    %check if memory permits all in one step for speed 
    if (length(master_indices)+length(slave_indices))*len*size(master_datarun.eis{master_indices(1)},1)<400*519*80
        
        if verbose, fprintf('\n map_ei...'); end
        
        mat = zeros(length(master_indices) + length(slave_indices), len * size(master_datarun.eis{master_indices(1)},1));
        for i=1:length(master_indices)
            t = master_datarun.eis{master_indices(i)}(:,1:len);
            if isempty(t), error('EI is empty'); end
            t(find(max(t,[],2)<electrode_threshold)) = 0; 
            mat(i,:) = t(:);
        end
        for i=1:length(slave_indices)
            t=slave_datarun.eis{slave_indices(i)}(:,1:len);
            if isempty(t), error('EI is empty');end
            t(find(max(t,[],2)<electrode_threshold))=0; 
            mat(length(master_indices)+i,:)=t(:);
        end

        %correlate
        temp=corrcoef(mat');
        corr=temp(1:length(master_indices),length(master_indices)+1:end);
        
        if verbose, fprintf('...finished\n'); end

    else
        
        if verbose, T=text_waitbar('map_ei'); end

        mat=zeros(length(master_indices)+1,len*size(master_datarun.eis{master_indices(1)},1));
        corr=zeros(length(master_indices),length(slave_indices));
        for i=1:length(master_indices) 
            t=master_datarun.eis{master_indices(i)}(:,1:len);
            if isempty(t), error('EI is empty');end
            t(find(max(t,[],2)<electrode_threshold))=0; 
            mat(i+1,:)=t(:);
        end 

        for i=1:length(slave_indices)
            t=slave_datarun.eis{slave_indices(i)}(:,1:len);
            if isempty(t), error('EI is empty');end
            t(find(max(t,[],2)<electrode_threshold))=0;
            mat(1,:)=t(:);

            temp=corrcoef(mat');
            corr(:,i)=temp(1,2:end);

            if verbose, T=text_waitbar(T,i/length(slave_indices)); end
        end
    end
end


%use only spatial info   
if space_only 
    exclude_list = [];

    if verbose, fprintf('\n map_ei...'); end
    mat = zeros(size(master_datarun.ei.eis{master_indices(1)},1),length(master_indices)+length(slave_indices));
    for i = 1:length(master_indices)
        st = max(master_datarun.ei.eis{master_indices(i)},[],2);
        temp = find(st < electrode_threshold);
        st(temp)=0;
        if length(find(st > electrode_threshold)) < significant_electrodes
            exclude_list = [exclude_list master_indices(i)];
        end
        mat(:,i) = st;
    end
    for i = 1:length(slave_indices)
        st = max(slave_datarun.ei.eis{slave_indices(i)},[],2);
        st(find(st < electrode_threshold)) = 0;
        mat(:,length(master_indices)+i) = st;
    end

    %correlate
    temp = corrcoef(mat);
    corr = temp(1:length(master_indices),length(master_indices)+1:end);

    if verbose, fprintf('...finished\n'); end
end


% If EI is blank, e.g. if this was a concatenated mapping run and the cell
% has no spikes in this run, then we may have EIs that are all zeros.
corr(isnan(corr)) = 0;


%sort
failed_cell_count = 0;
cell_list_map = cell(1,length(master_datarun.cell_ids));
for i = 1:length(master_indices)
    %[largest_corr, max_corr_index] = max(corr(i,:));    
    [t_corr, t_corr_index] = sort(corr(i,:));
    largest_corr=t_corr(end);
    max_corr_index=t_corr_index(end);
    
    if largest_corr >= corr_threshold 
        [largest_corr_col, max_corr_index_col] = max(corr(:,max_corr_index));
        if largest_corr==largest_corr_col
            cell_list_map{master_indices(i)} = slave_datarun.cell_ids(slave_indices(max_corr_index));
            if troubleshoot
                fprintf('master cell_ID:%4d -> %4d (correlation:%1.3f) - next closest: %4d(%1.3f) %4d(%1.3f) %4d(%1.3f) - distance: %2.2fsd\n', master_datarun.cell_ids(master_indices(i)), cell_list_map{master_indices(i)}, largest_corr,...
                    slave_datarun.cell_ids(slave_indices(t_corr_index(end-1))),t_corr(end-1), slave_datarun.cell_ids(slave_indices(t_corr_index(end-2))),t_corr(end-2), slave_datarun.cell_ids(slave_indices(t_corr_index(end-3))),t_corr(end-3),(largest_corr-t_corr(end-1))/robust_std(corr(i,:)) )
                %pause(0.1)
            end
        else
        	fprintf('master cell_ID:%4d best slave %4d (correlation:%1.3f) taken by master %4d(%1.3f)\n', master_datarun.cell_ids(master_indices(i)), slave_datarun.cell_ids(slave_indices(max_corr_index)), largest_corr,...
                    master_datarun.cell_ids(master_indices(max_corr_index_col)),largest_corr_col )
                output_matrix(i,:) = [master_datarun.cell_ids(master_indices(i)), slave_datarun.cell_ids(slave_indices(max_corr_index)), master_datarun.cell_ids(master_indices(max_corr_index_col))];
        end
    else
        if troubleshoot
            fprintf('no match for cell %d ', master_datarun.cell_ids(master_indices(i)))
            fprintf('- the highest correlation was %1.3f for cell %d \n', largest_corr,slave_datarun.cell_ids(slave_indices(max_corr_index))) 
            %pause(0.1)
        end
        failed_cell_count = failed_cell_count + 1;
        failed_cells(failed_cell_count) = master_datarun.cell_ids(master_indices(i));
        cell_list_map{master_indices(i)}=[];
    end
end















