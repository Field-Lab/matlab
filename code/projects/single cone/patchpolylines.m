function p = patchpolylines(lines, varargin)
% PATCHPOLYLINES
% usage: p = patchpolylines(lines, opts)
%
% Plots output from SEGS2POLY where SEGS is, for example, output from MAP2MANHATTAN
%
% 2012-06 phli
%

opts = inputParser();
opts.addParamValue('colors', []);
opts.addParamValue('cmf', @jet)
opts.addParamValue('scale', [1 1]);
opts.addParamValue('patchopts', {});
opts.addParamValue('fill', false);
opts.addParamValue('fillcolors', []);
opts.parse(varargin{:});
opts = opts.Results;

n = length(lines);

if isempty(opts.colors),                  opts.colors = opts.cmf(n); end
if opts.fill && isempty(opts.fillcolors), opts.fillcolors = opts.colors; end
if ~isempty(opts.fillcolors),             opts.fill = true; end


oldhold = ishold();
if ~oldhold, cla; end
hold on;

for i = 1:n
    color = opts.colors(i,:);
    if opts.fill, facecolor = opts.fillcolors(i,:);
    else facecolor = 'none'; end
    
    linesi = lines{i};
    nlines = length(linesi);
    
    % Setup for multiobject patch plot
    longest = max(cellfun(@(l)(size(l,2)), linesi));
    xdata = NaN(longest, nlines);
    ydata = xdata;    
    
    % Load data into multiobject patch array; shorter rows will remain
    % padded with NaNs at end
    for j = 1:nlines
        line = linesi{j};
        xdata(1:length(line),j) = line(1,:)';
        ydata(1:length(line),j) = line(2,:)';
    end

    % Rescale and patch
    xdata = (xdata-0.5) .* opts.scale(1) + 0.5;
    ydata = (ydata-0.5) .* opts.scale(2) + 0.5;
    p(i) = patch(xdata, ydata, color, 'EdgeColor', color, 'FaceColor', facecolor, opts.patchopts{:});
end
axis equal tight;

if ~oldhold, hold off; end
if nargout < 1, clear p; end