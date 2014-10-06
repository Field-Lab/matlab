function varargout = find_local_maxima(matrix, varargin)
% find_local_maxima     identify entries which are the largest in their neighborhood
%
% usage:  maxes = find_local_maxima(matrix);
%         indices = find_local_maxima(matrix,'return','indices');
%         coordinates = find_local_maxima(matrix,'return','coordinates');
% 
%
% arguments:   matrix - double matrix
%            varargin - struct or list of optional parameters (see below)
%
% outputs:     depends on value of 'return'
%       'matrix'  - logical matrix showing local maxima
%       'indices' - matrix indices of local maxima  (row,column)
%   'coordinates' - coordinates of local maxima  (x,y)
%
%
%
% optional params, their default values, and what they specify:
%
% radius           	1               radius of the region in which an entry must be the local maximum
% thresh            0               only search entries greater than this threshold (reduces computation time)
% return            'matrix'        see outputs above
% verbose           false           show progress
%
% 2008-10
% 2012-07 phli, Check for off edge needed to be changed; intersect very slow for large matrices
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('radius',1,@(x) (mod(x,1)==0) && (x > 0)); % positive integer
p.addParamValue('thresh',0);
p.addParamValue('return', 'matrix', @(x)any(strcmpi(x,{'matrix','indices','coordinates'})));
p.addParamValue('verbose',false);

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;




% BODY OF THE FUNCTION



% radius of search
search_rad = params.radius;

% initialize matrix to store significant stixels
local_maxima = false([size(matrix,1) size(matrix,2)]);

% only search pixels that exceed threshold
search_locations = matrix > params.thresh;

% don't search around edges
%search_locations([1:search_rad end-search_rad+1:end],:) = false;
%search_locations(:,[1:search_rad end-search_rad+1:end]) = false;

% get coordinates of search locations
[search_y search_x] = find(search_locations);
 
% show progress
if params.verbose
    T=text_waitbar('Searching for local maxima...');
end

% go through each one, saving the locations of local maxima
for cc = 1:length(search_y)

    %temp=zeros(size(search_locations));temp(search_y(cc),search_x(cc))=1;
    %figure(60);clf;imagesc(temp)

    % note search range
    y_range = search_y(cc) + [-search_rad:search_rad];
    x_range = search_x(cc) + [-search_rad:search_rad];
       
    % be sure not to try to search off the edge
    y_range = y_range(y_range > 1 & y_range < size(matrix,1));
    x_range = x_range(x_range > 1 & x_range < size(matrix,2));

    % get local neighborhood
    local_vals = matrix(y_range,x_range);
    
    % note location of original center
    orig_cent_y = find(search_y(cc)==y_range);
    orig_cent_x = find(search_x(cc)==x_range);

    % locate pixels that are the maximum value
    max_vals = (local_vals == max(max(local_vals)));

    switch 'center'
        case 'any' % save local maximum, wherever it's located
            % save these values in the larger local_maxima
            local_maxima(y_range,x_range) = local_maxima(y_range,x_range) + max_vals;

        case 'center' % save it only if it is in the center
            if max_vals(orig_cent_y,orig_cent_x)
                local_maxima(search_y(cc),search_x(cc)) = true;
            end
    end
    
    % show progress
    if params.verbose && mod(cc,1000)==0
        T=text_waitbar(T,cc/search_y);
    end
   
end


switch params.return
    case 'matrix'
        varargout = {local_maxima};
    case 'indices'
        [i,j] = find(local_maxima);
        varargout = {[i j]};
    case 'coordinates'
        [i,j] = find(local_maxima);
        varargout = {[j i]};
end


