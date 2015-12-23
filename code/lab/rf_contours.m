function [contour_polygons,extras,params] = rf_contours(summary, contour_levels, params)
% rf_contours     compute contours for an sta spatial summary
%
% usage:  [contour_polygons,extras,params] = rf_contours(summary, contour_levels, params)
%
% arguments:        summary - STA spatial summary, 2-d or 3-d matrix (height, width, color)
%            contour_levels - vector giving contour levels
%                    params - struct of optional parameters (see below)
%
% outputs:   contour_polygons - cell array of contours at each specified level
%                                   contours are in polygon format (see PolygonClipper)
%                      extras - struct with various intermediate results (see below)
%                      params - struct with parameters of this function (see below)
%
%
% fields in extras and what they contain:
%
%   extras.outer_bounds     outer boundaries of contours at each contourlevel
%
%
% optional fields in params, their default values, and what they specify:
%
% color       	[]  	which color to use.  if blank, use color with highest peak.
%
%


% SET UP OPTIONAL ARGUMENTS

% if not specified, make params empty
if ~exist('params','var');params = [];end

% specify default parameters
defaults.color = [];

% combine user and default parameters
params = default_params( defaults, params);



% DEAL WITH SPECIAL CASE

% if empty or constant, return empties
if isempty(summary) || max(max(max(summary))) == min(min(min(summary)))
    contour_polygons = [];
    extras.outer_bounds = [];
    params = params;
    return
end



% REDUCE TO ONE COLOR

% set color to 1 if there's only one color channel
if size(summary,3) == 1
    params.color = 1;
end

% otherwise, if a color was not specified, pick the one with the greatest peak
if isempty(params.color)
    % for each color, get the greatest deviation from 0
    color_peaks = squeeze(max(max(max(summary,[],1),[],2),[],4));
    % get the maximum deviation, and find which colors attain this value
    max_colors = find(color_peaks == max(color_peaks));
    % get the index of this color
    % if there's a tie (possible but very unlikely), just use the first color
    params.color = max_colors(1);
end

% save only the desired color
summary = summary(:,:,params.color);


% COMPUTE CONTOURS

% go through each contour level
for cc = 1:length(contour_levels)

    % compute contours in stupid Matlab format
    contour_temp = contourc(summary,contour_levels(cc) * [1 1]);

    % convert to polygon format and save to cell array
    [contour_polygons{cc},extras.outer_bounds{cc}] = to_polygon_format(contour_temp);
    
    % Record level info directly in struct
    [contour_polygons{cc}.level] = deal(contour_levels(cc));
end





function [contour_polygons,outer_bounds] = to_polygon_format(cnt)
% put each cell's contours into the polygon struct format


% deal with empty case
if isempty(cnt)
    % return an empty polygon
    contour_polygons = struct('x',[0],'y',[0],'hole',0);
    outer_bounds.x = [0 0];
    outer_bounds.y = [0 0];
    return
end



% TRANSFORM THE MATLAB CONTOURS TO A MORE REASONABLE FORMAT

% find the number of contour bands (where a "contour band" is a separate polygon)
band_start_indices = find( cnt(1,:) == cnt(1,1) );
num_bands = length(band_start_indices);

% load the bands into a cell array: contour_bands{1}, contour_bands{2},...
band_start_indices = [band_start_indices size(cnt,2)+1];
for bb = 1:num_bands
    start = band_start_indices(bb) + 1;
    finish = band_start_indices(bb+1) - 1;
    contour_bands{bb} = cnt(1:2,start:finish);
end



% GET THE OUTER BOUNDS OF THESE CONTOURS

% get all points
contour_points = [];
for bb = 1:length(contour_bands)
    contour_points = [contour_points contour_bands{bb}];
end
% identify their extreme points
outer_bounds.x = [min(contour_points(1,:)) max(contour_points(1,:))];
outer_bounds.y = [min(contour_points(2,:)) max(contour_points(2,:))];



% IDENTIFY WHICH CONTOUR BANDS CONTAIN WHICH OTHERS

% if there's just contour band, this is easy
if length(contour_bands) == 1
    contour_polygons = struct('x',contour_bands{1}(1,:),'y',contour_bands{1}(2,:),'hole',0);

else
    % otherwise, do work :(

    % go through each contour band
    for bb = 1:length(contour_bands)

        % determine how many other bands contain this one
        num_containing_polygons = 0;
        for pp = 1:length(contour_bands)
            num_containing_polygons = num_containing_polygons + ...
                polygon_contained_in_polygon(contour_bands{bb},contour_bands{pp},[1 0 0]);
        end
        % account for the fact that the polygon contains itself
        num_containing_polygons = num_containing_polygons - 1;

        % assign the polygon sign
        if mod(num_containing_polygons,2)==0
            % encloses positive area
            hole_value = 0;
        else
            % encloses negative area
            hole_value = 1;
        end

        % load the polygon
        contour_polygons(bb) = struct('x',contour_bands{bb}(1,:),'y',contour_bands{bb}(2,:),'hole',hole_value);

    end
end




function contains = polygon_contained_in_polygon(pg1,pg2,return_values)
% pg1 is contained in pg2 returns return_values(1)
% pg2 is contained in pg1 returns return_values(2)
% neither returns return_values(3)
% pg2 == pg1 returns return_values(1)

pg1INpg2 = inpolygon(pg1(1,:),pg1(2,:),pg2(1,:),pg2(2,:));
pg2INpg1 = inpolygon(pg2(1,:),pg2(2,:),pg1(1,:),pg1(2,:));

if sum(pg1INpg2) == length(pg1INpg2)
    contains = return_values(1);
else
    if sum(pg2INpg1) == length(pg2INpg1)
        contains = return_values(2);
    else
        contains = return_values(3);
    end
end




