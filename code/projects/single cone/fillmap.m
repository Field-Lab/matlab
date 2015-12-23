function h = fillmap(map, varargin)
% FILLMAP   Create filled squares for each pixel in the given bitmap, colors by index
%
% 2012-06 phli
%

opts = inputParser();
opts.addParamValue('colors', []);
opts.addParamValue('cmf', @jet)
opts.addParamValue('scale', [1 1]);
opts.addParamValue('fillopts', {});
opts.parse(varargin{:});
opts = opts.Results;


indices = setdiff(unique(map), 0);
n = length(indices);

if isempty(opts.colors)
    opts.colors = opts.cmf(n);
end


fills = cell(1,3*n);
for i = 1:n
    index = indices(i);
    color = opts.colors(i,:);

    [r c] = find(map == index);
    x = [c'-0.5; c'+0.5; c'+0.5; c'-0.5];
    y = [r'-0.5; r'-0.5; r'+0.5; r'+0.5];
    
    x = (x-0.5).*opts.scale(1) + 0.5;
    y = (y-0.5).*opts.scale(2) + 0.5;
    
    fills = [fills x y color];
end
h = fill(fills{:}, 'EdgeColor', 'none', opts.fillopts{:});
axis equal tight;

if nargout < 1, clear h; end