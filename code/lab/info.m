function info(datarun, params)
% INFO     displays overview of datarun
%
% usage:  info(datarun, params)
%
% arguments:  datarun - datarun struct with field specifying X
%              params - display_cell_ids, number of displayed cell_ids,
%              default 10, -1 displays all 
%                       
% greschner


% SET UP OPTIONAL ARGUMENTS
% if not specified, make params empty
if ~exist('params','var');params = [];end
% specify default parameters
defaults.display_cell_ids = 10;
% combine user and default parameters
params = default_params( defaults, params);


for j=1:length(datarun)
    if iscell(datarun)
        data=datarun{j};
    else
        data=datarun;
    end
      
    if isfield(data.names,'short_name')
        fprintf('\ndatarun{%d}: %s\n',j,data.names.short_name);
    end
        
    if isfield(data,'cell_types')
        for i=1:length(data.cell_types)
            fprintf('\n%3.0f.  (%3.0f)  %s:      ',i,length(data.cell_types{i}.cell_ids),data.cell_types{i}.name);   

            temp=length(data.cell_types{i}.cell_ids);
            if params.display_cell_ids==-1 | params.display_cell_ids>temp
                display_cell_ids=temp;
            else
                display_cell_ids=params.display_cell_ids;
            end

            for ii=1:display_cell_ids
                fprintf(' %d',data.cell_types{i}.cell_ids(ii));
            end

        end
    end
    fprintf('\n');
end

