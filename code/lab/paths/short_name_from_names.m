function short_name = short_name_from_names(names)
% short_name_from_names     generate short name for a dataset from the short name
%
% usage:  short_name = short_name_from_names(names)
%
% arguments:     names - names struct from datarun.names
%
% outputs:     short_name -
%
%
%
% gauthier 2008-10
%







% BODY OF THE FUNCTION


% use experiment and condition if possible
if isfield(names,'experiment') && isfield(names,'condition')
    short_name = [names.experiment '_' names.condition];

else
    % cycle through various options
    field_names = {'rrs_prefix','rrs_params_path','rrs_neurons_path','rrs_ei_path','rrs_sta_path'};
    
    for ff = 1:length(field_names)
        if isfield(names,field_names{ff}) && ~isempty(names.(field_names{ff}))
            short_name = strrep(names.(field_names{ff}),server_path,'');
            break
        end
    end
    
    % remove file extension, if present
    if exist('short_name','var') && ff > 1
        % get locations of all periods
        period_index = strfind(short_name,'.');
        % remove everything after the last one
        short_name = short_name(1:period_index(end)-1);
    end
end

% if a name has been given...
if exist('short_name','var')
    % clean it up
    short_name = strrep(short_name,'/','_');
else
    % otherwise, assign a random name
    short_name = sprintf('unknown_dataset_%d',round(sum(clock)));
end
