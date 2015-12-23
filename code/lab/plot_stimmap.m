function plot_stimmap(map, varargin)

opts = inputParser();
opts.addParamValue('EdgeColor', []);
opts.addParamValue('FaceColor', 'none');
opts.addParamValue('collapse_indices', false);
opts.addParamValue('plot_image', true);
opts.parse(varargin{:});
opts = opts.Results;

if isempty(opts.EdgeColor)
    if opts.plot_image
        opts.EdgeColor = 'w';
    else
        opts.EdgeColor = 'k';
    end
end


oldhold = ishold();
if ~oldhold, cla; end
hold on;


if opts.plot_image
    imagesc(map);
end


[edges_x, edges_y] = map2manhattan(map, 'collapse_output', opts.collapse_indices);
for i = 1:length(edges_x)
    if iscell(opts.EdgeColor), edgecolor = opts.EdgeColor{i}; 
    else                       edgecolor = opts.EdgeColor; end
    
    if iscell(opts.FaceColor), facecolor = opts.FaceColor{i};
    else                       facecolor = opts.FaceColor; end
    
    patch(edges_x{i}, edges_y{i}, edgecolor, 'FaceColor', facecolor, 'EdgeColor', edgecolor);
end

axis equal;
axis([0 size(map,1) 0 size(map,2)]);


if ~oldhold, hold off; end