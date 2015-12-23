function params = default_params( defaults, params)
% DEFAULT_PARAMS     combine default parameters with user specified parameters
%
% usage:  params = default_params( defaults, params);
%
% arguments:    defaults - struct of default parameters
%                 params - struct of parameters supplied by user
%
% outputs:        params - struct combining the default values with user values
%
% Returns an error if any field in params is not found in defaults, and prints
% the list of valid fields.
%
% gauthier, 2008-03-01
% greschner ...  || isnan(params.(default_param_list{pp}))


% if params is empty, just set it to defaults
if isempty(params)
    params = defaults;
    return
end


% get list of default parameters
default_param_list = fieldnames(defaults);

% get list of user parameters
user_param_list = fieldnames(params);


% check if any fields of params are not found in defaults
for pp = 1:length(user_param_list)
    if ~isfield(defaults,user_param_list{pp})
        % give an error and show the list of valid arguments
        
        %error('Optional argument ''%s'' not recognized.  Valid arguments are %s.',...
        %    user_param_list{pp},list_from_cell(default_param_list))
        fprintf('\noptional parameters and default values:\n\n')
        disp(defaults)
        error('Optional argument ''%s'' not recognized.  Valid arguments and default values are printed above.',user_param_list{pp})
    end
end


% cycle through default parameters...
for pp = 1:length(default_param_list)
    % keep any user entered value...
    if ~isfield(params,default_param_list{pp}) % || isequal('default value',params.(default_param_list{pp}))
        % but otherwise use the default value
        params.(default_param_list{pp}) = defaults.(default_param_list{pp});
    end
end



function list = list_from_cell(fn)
% takes in field names, returns single string listing the names separated by commas

list = ['''' fn{1} ''''];
for ff = 2:(length(fn)-1)
    list = [list ', ''' fn{ff} ''''];
end
list = [list ', and ''' fn{end} ''''];
