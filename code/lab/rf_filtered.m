function [filt_summary, params] = rf_filtered(summary, params)
% rf_filtered     Filter an STA spatial summary
%
% usage:  result = rf_filtered(summary, params)
%
% arguments:  summary - STA spatial summary, 2-d or 3-d matrix (height, width, color) 
%              params - struct of optional parameters (see below)
%
% outputs:    filt_summary - filtered version of spatial summary
%
%
% optional fields in params, their default values, and what they specify:
%
% filt_type         'gauss'         type of filter to apply.  see options below.
%                                       if empty, reverts to default
% filter_size    	[]              sidelength of filter matrix
%                                       if empty, determined by the specific filter parameters
%
%
% several options in params are specific to a particular kind of filter:
%
%   'gauss'     perform gaussian blur
%
%       radius     	3       radius of the blur
%
%
%   'none'      perform no filtering
%
%
%
%
% % greschner 'none'



% SET UP OPTIONAL ARGUMENTS

% if not specified, make params empty
if ~exist('params','var');params = [];end

% specify default parameters
defaults.filt_type = 'gauss';
defaults.filter_size = [];
defaults.radius = 3;

% combine user and default parameters
params = default_params( defaults, params);

% if filter type field is empty, revert to default
if isempty(params.filt_type)
    params.filt_type = defaults.filt_type;
end


% hand special case where summary is empty
if isempty(summary)
    filt_summary = [];
end


% APPLY FILTER

switch params.filt_type
    case 'gauss'
        
        % get filter parameters
        gauss_radius = params.radius;
        if isempty(params.filter_size)
            filter_size = max(5,round(3*gauss_radius));
        else
            filter_size = params.filter_size;
        end
        
        % compute filter
        summ_filt = fspecial('gaussian',filter_size,gauss_radius);
        
        % apply filter
        filt_summary = imfilter(summary,summ_filt/sum(sum(summ_filt)));
    
    case 'none' 
        disp('rf_filtered: not filtered') 
        filt_summary = summary;
    otherwise
        error('rf_filtered: Filter type ''%s'' not recognized.',params.filt_type)
end

