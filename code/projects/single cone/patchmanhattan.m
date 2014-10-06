function p = patchmanhattan(nyc, varargin)
% PATCHMANHATTAN
% usage: p = patchmanhattan(nyc, opts)
%
% Plots output from MAP2MANHATTAN
%
% 2012-02 phli
%

opts = inputParser();
opts.addParamValue('colors', []);
opts.addParamValue('cmf', @jet)
opts.addParamValue('scale', [1 1]);
opts.addParamValue('patchopts', {});
opts.parse(varargin{:});
opts = opts.Results;

if ~iscell(nyc), nyc = {nyc}; end
n = length(nyc);

if isempty(opts.colors)
    opts.colors = opts.cmf(n);
end

oldhold = ishold();
if ~oldhold, cla; end
hold on;

for i = 1:n
    colormap = opts.colors(i,:);
    p(i) = patch((nyc{i}.x - 0.5) .* opts.scale(1) + 0.5, (nyc{i}.y - 0.5) .* opts.scale(2) + 0.5, colormap, 'EdgeColor', colormap, opts.patchopts{:});
end
axis equal tight;

if ~oldhold, hold off; end
if nargout < 1, clear p; end