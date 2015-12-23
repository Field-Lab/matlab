function [ccf time]=get_correlation_shift(datarun, cell_specification_1, cell_specification_2, varargin) 


% SET UP OPTIONAL ARGUMENTS
    p = inputParser;
    p.addParamValue('dt', .004);
    p.addParamValue('offset', .15);
    p.addParamValue('shift', []);
    p.addParamValue('stop', []);
    p.addParamValue('start', []);
    p.addParamValue('max_isi', 0);
    p.parse(varargin{:});
    params = p.Results;
         
       
 
% get cell numbers
index1=get_cell_indices(datarun,cell_specification_1);
index2=get_cell_indices(datarun,cell_specification_2);

spikes_1=datarun.spikes{index1};
spikes_2=datarun.spikes{index2};

[ccf time]=correlation_shift(spikes_1, spikes_2, params);


    
    
    






























    
    