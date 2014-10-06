function [com, extras, params] = rf_com(sta, varargin)
% rf_com     Get the center of mass (COM) of an RF frame
%
% usage:  [com, extras, params] = rf_com(sta, varargin)
%                               
% arguments:       sta - YxXxZxC double matrix (C = # colors) representing the RF
%                  rf = a YxXxC rf can be handed to rf_com instead of an sta
%            varargin - struct or list of optional parameters (see below)
%
% outputs:        com - coordinates of the COM in the form [x y]
%              extras - intermediate results of the computation
%              params - parameters used in the calculation
%
%
% optional fields in params, their default values, and what they specify:
%
% sig_stixels       []          significant stixels to use.  If emtpy, they are computed from scratch.
% ss_params         struct     	struct of parameters to pass to function significant_stixels (if sig_stixels is empty)
%                                   see significant_stixels for more information
%
% exclude_outliers  true        stixels that are significant but have a mean distance > 'distance_thresh' from other
%                               significant stixels are excluded
% distance_thresh   50          the threshold distance for a sig. stixel being defined as an outlier
%
% positive          false       if true, then sig_stixs must have the polarity of the strongest stixel
%
%
% output fields in extras:
%
%             none yet!
%
% 2008-10 gauthier
%    edited: 2009-08  gdf (added excluded_outliers and distance_thresh)
%            2009-08  gdf (computation now takes stixel weights into account)
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('ss_params', struct);
p.addParamValue('sig_stixels', []);
p.addParamValue('distance_thresh', 50, @isnumeric);
p.addParamValue('exclude_outliers', false, @islogical);
p.addParamValue('positive', false, @islogical)

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;



% check sta has appropriate dimesnions for sta or rf
array_size = length(size(sta));
if array_size == 4
    rf = rf_from_sta(sta);
elseif (array_size == 3) || (array_size == 2)
    rf = sta;
else
    error('Check Inputs: sta has neither the dimensions of an sta nor an rf.  size of sta is [%s]',num2str(size(sta)))
end


% IDENTIFY SIGNIFICANT STIXELS

% if user did not provide significant stixels
if isempty(params.sig_stixels)
    % get signficant stixels
    if array_size == 4
        [sig_stixels, ss_params] = significant_stixels(sta,params.ss_params);
    else
        [sig_stixels, ss_params] = significant_stixels(rf,params.ss_params);
    end
    
else
    % if user did provide them, use them
    sig_stixels = params.sig_stixels;
    %set ss_params to empty struct
    ss_params = struct;
end

% warn if there are no significant stixels
temp_indices = find(full(sig_stixels));
if isempty(temp_indices)
    warning('no significant stixels found')
end  

% exclude outlier stixels
if params.exclude_outliers
    [temp_rows, temp_cols] = ind2sub(size(sig_stixels), temp_indices);

    % get locations of stixels and calculate their distances
    locations = [temp_rows, temp_cols];
    distances = ipdm(locations);
    mean_distances = mean(distances, 1);

    % exclude stixels that are further away than threshold
    exclude_indices = find(mean_distances > params.distance_thresh);

    sig_stixels(temp_rows(exclude_indices), temp_cols(exclude_indices)) = false;
end

% ensure values at sig_stixels are postive
if params.positive
    weights = squeeze(sum(rf,3)) .* sig_stixels;

    keeper_indices = find(weights > 0);
    [temp_rows, temp_cols] = ind2sub(size(weights), keeper_indices);

    pos_stixels = false(size(weights));
    pos_stixels(temp_rows, temp_cols) = true;

    sig_stixels = sig_stixels & pos_stixels;
    sig_stixels = sparse(sig_stixels);
end    

% handle case where there are no significant stixels
if isempty(sig_stixels) || max(max(sig_stixels)) == 0
    com = [];
    extras = [];
    params = struct;
    return
end

% note parameters
params = ss_params;


% COMPUTE COM USING ONLY THE SIGNIFICANT STIXELS

% compute the COM of stixels greater than thresh
gs_rf = squeeze(sum(rf, 3));
% only keep positive pixels (negative ones throw off the calculation)
gs_rf = gs_rf .* (gs_rf > 0);
com = getfield(regionprops(double(full(sig_stixels)), gs_rf,'WeightedCentroid'),'WeightedCentroid');


% return other information if requested
if nargout > 1
    extras = [];
end

