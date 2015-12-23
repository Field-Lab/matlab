function plot_voronoi(datarun, varargin)
% 2011-07 phli

opts = inputParser;
opts.addParamValue('mode', 'polygons');
opts.addParamValue('bounds', []);
opts.addParamValue('plot_centers', false);
opts.addParamValue('label_cones', false);
opts.addParamValue('label_color', 'y');
opts.addParamValue('border_color', 'k');
opts.addParamValue('center_color', 'k');
opts.addParamValue('highlight_regions', []);
opts.addParamValue('highlight_rgba', []);
opts.addParamValue('highlight_legend', {});
opts.addParamValue('highlight_ranked_cones', []);
opts.addParamValue('cmf', @jet);
opts.addParamValue('cone_rank_rgc', []);
opts.addParamValue('mask_space_neighbor_min', []);
opts.addParamValue('mask_space_self_max', []);
opts.parse(varargin{:});
opts = opts.Results;


oldhold = ishold;
hold on;


% Highlight selected cones; may select directly by index or by weight rank
if ~isempty(opts.highlight_ranked_cones)
    ranked_cones = get_cones_by_weight_rank(datarun, opts.cone_rank_rgc, opts.highlight_ranked_cones);
    opts.highlight_regions = [opts.highlight_regions ranked_cones];
end

if ~isempty(opts.highlight_regions) && isempty(opts.highlight_rgba)
    opts.highlight_rgba = opts.cmf(length(opts.highlight_regions));
end


centers = datarun.cones.centers;
if ~isempty(opts.bounds)
    inxbounds = centers(:,1) > opts.bounds(1) & centers(:,1) < opts.bounds(2);
    inybounds = centers(:,2) > opts.bounds(3) & centers(:,2) < opts.bounds(4);
    inbounds = inxbounds & inybounds;
    
    centers = centers(inbounds, :);
    
    reindexed_regions = [];
    for i = 1:numel(opts.highlight_regions)
        oldindex = opts.highlight_regions(i);
        newindex = find(find(inbounds) == oldindex, 1);
        if ~isempty(newindex);
            reindexed_regions(end+1) = newindex;
        end
    end
    opts.highlight_regions = reindexed_regions;
end

switch opts.mode
    case 'polygons'
        highlights = plot_voronoi_polygons(centers, opts);

    case 'manhattan'
        if ~isempty(opts.mask_space_neighbor_min) || ~isempty(opts.mask_space_self_max)
            masks = space_voronoi_masks(datarun, opts.mask_space_neighbor_min, opts.mask_space_self_max);
        else
            masks = datarun.cones.mosaic.voronoi_masks;
        end
        if ~isempty(opts.bounds), masks = masks(inbounds); end
        scale = 1 / datarun.cones.mosaic.voronoi_masks_scale;
        highlights = plot_voronoi_manhattan(masks, scale, opts);
end

if ~isempty(opts.highlight_legend)
    legend(highlights, opts.highlight_legend{:});
end


if opts.label_cones
    for i = 1:size(datarun.cones.centers, 1)
        if ~isempty(opts.bounds) && ~inbounds(i), continue; end
        t = text(datarun.cones.centers(i,1), datarun.cones.centers(i,2), num2str(i));
        set(t, 'Color', opts.label_color);
    end
end


if opts.plot_centers
    h = plot(centers(:,1), centers(:,2), '.');
    set(h, 'color', opts.center_color);
end

axis equal;
if ~isempty(opts.bounds), axis(opts.bounds); end

if ~oldhold, hold off; end



function highlights = plot_voronoi_manhattan(masks, scale, opts)

disp(['Calculating ' num2str(length(masks)) ' Manhattan borders.']);
for i = 1:length(masks)
    fprintf('=');
end
fprintf('\n');

p = zeros(length(masks), 1);
for i = 1:length(masks)
    fprintf('*');
    [edges_x,edges_y] = mask2manhattan(masks{i});
    % gdf: removed this edge correction which is thought to cause errors
    %edges_x = ((edges_x - 0.5) .* scale) + 0.5;
    %edges_y = ((edges_y - 0.5) .* scale) + 0.5;
    edges_x = edges_x .* scale;
    edges_y = edges_y .* scale;
    if isempty(edges_x), continue; end;
    p(i) = patch(edges_x, edges_y, opts.border_color, 'FaceColor', 'none', 'EdgeColor', opts.border_color);
end
fprintf('\n');

highlights = [];
for i = 1:length(opts.highlight_regions)
    regionh = p(opts.highlight_regions(i));
    rgb = opts.highlight_rgba(i,1:3);

    if size(opts.highlight_rgba, 2) == 4
        alpha = opts.highlight_rgba(i,4);
        set(regionh, 'FaceColor', rgb, 'FaceAlpha', alpha);
    else
        % Just highlight border
        set(regionh, 'EdgeColor', rgb, 'LineWidth', 2);
    end
    
    uistack(regionh, 'top');
end



function highlights = plot_voronoi_polygons(centers, opts)
highlights = [];

% Plot tesselation
if ~strcmp(opts.border_color, 'none')
    [vx,vy] = voronoi(centers(:,1), centers(:,2));
    lines = plot(vx, vy, opts.border_color);
    set(lines, 'xliminclude', 'off', 'yliminclude', 'off')
end
if isempty(opts.highlight_regions), return; end

% Highlight specified regions
hold on
[v,c] = voronoin(centers);

for i = 1:length(opts.highlight_regions)
    region = opts.highlight_regions(i);
    rgb   = opts.highlight_rgba(i,1:3);

    if size(opts.highlight_rgba, 2) == 4
        alpha = opts.highlight_rgba(i,4);
        highlights(i) = patch(v(c{region}, 1), v(c{region}, 2), rgb, 'FaceAlpha', alpha);
    else
        % Just highlight border
        highlights(i) = patch(v(c{region}, 1), v(c{region}, 2), rgb, 'FaceColor', 'none', 'EdgeColor', rgb, 'LineWidth', 2);
    end
end