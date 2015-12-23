function datarun=correlation(datarun, cell_specification_1, cell_specification_2, varargin) 
%returns matrix of cell_specification_1 x cell_specification_2 with the
%correlation strength log2(p_ab/(p_a*p_b))
%
% synchrony_index=correlation(datarun, cell_specification_1, [cell_specification_2], [params]) 
%
% defaults.bin = 10; 


% SET UP OPTIONAL ARGUMENTS
    p = inputParser;
    p.addParamValue('synchrony_index_field', 'synchrony_index');% 
    p.addParamValue('corrcoef', 1);%
    p.addParamValue('bin', 10);%
    p.addParamValue('shift', 0);%
   
    % parse inputs
    p.parse(varargin{:});
    params = p.Results;
  
        
%get synchrony_index        
if isfield(datarun, params.synchrony_index_field);
    synchrony_index=datarun.(params.synchrony_index_field);
else
    synchrony_index=zeros(length(datarun.cell_ids));
    synchrony_index=sparse(synchrony_index);
end
        
        
% get cell numbers
index_1=get_cell_indices(datarun,cell_specification_1);
if ~exist('cell_specification_2','var');
    index_2=index_1;
else
    index_2 = get_cell_indices(datarun,cell_specification_2);
end



if isequal(index_1, index_2)
    sync=correlation_(datarun.spikes(index_1), [],'corrcoef',params.corrcoef,'bin',params.bin); 
else
    sync=correlation_(datarun.spikes(index_1), datarun.spikes(index_2),'corrcoef',params.corrcoef,'bin',params.bin,'shift',params.shift); 
end

synchrony_index(index_1, index_2)=sync;
synchrony_index(index_2, index_1)=sync';


datarun.synchrony_index=synchrony_index; 





if 0%qick test

    clf

    cell_type_1=1;
    cell_type_2=1;
    
    datarun{2}=correlation(datarun{2},{cell_type_1},{cell_type_2});
    [dist_stixel dist_sd normdist] = get_sta_fit_distance(datarun{1}, {cell_type_1},{cell_type_2});

    index_1=get_cell_indices(datarun{2},{cell_type_1});
    index_2=get_cell_indices(datarun{2},{cell_type_2});
    syn=full(datarun{2}.synchrony_index(index_1,index_2));

    plot(normdist,syn,'k.');

    set(gca,'XLim',[0 20]);    

end


