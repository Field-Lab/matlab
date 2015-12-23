function varargout = flexsubplot(traces, positions, varargin)
% FLEXSUBPLOT    Plots multiple traces custom positioned on a single set of axes
% usage: fig = flexsubplot(traces, positions, options)
%
% inputs: TRACES      Cell array of traces to plot.
%
%         POSITIONS   Nx2 matrix of x/y coordinates for centering location 
%                     of each plot.
%
% options:
%         FIG         gcf       Used for plotting
%
%         WIDTH       []        How wide each sub plot should be, in same
%                               units as POSITIONS.  If left blank, gets
%                               calculated automatically from the minimum
%                               x separation of POSITIONS.
%
%         HEIGHT      []        How high sub plots should be.  If YMODE is
%                               'local', every subplot ends up this height.
%                               If YMODE is 'global', then this is the
%                               maximum height for the plot that has the
%                               greatest amplitudes.  Often, it is useful
%                               to set this to something greater than the
%                               actual spacing of the subplots, so that a
%                               few of the traces will actually extend
%                               beyond their own window but most will not.
%                               Same units as POSITIONS.  If left blank,
%                               gets calculated automatically from the
%                               minimum y separation of POSITIONS.
%
%         YSCALE      1         Lets HEIGHT be calculated automatically but
%                               then scale it.  Most useful for scaling up
%                               so that the largest traces extend beyond
%                               their y windows.
%
%         YMODE       'global'  Can be set to 'global' or 'local'.  Global
%                               means the y axes are scaled the same for
%                               all traces.  Local means each y axis is
%                               scaled individually to fit its own subplot
%                               window.
%
%         TRACESX     []        Cell array of x-coordinates for traces.  Can be left
%                               blank for a given trace, or can be left entirely
%                               blank.
%
%         PLOTARGS    {{}}      Cell array of N cell arrays containing extra 
%                               arguments passed along to PLOT().  If some are
%                               left out, they will take their args from the
%                               previous index.
%
%         COLORS      []        Series of colors to cycle through.
%
%         COLOR_PROP  'Color'   The property to set to the cycle color on
%                               plot command; defaults to 'Color', which is 
%                               appropriate for 'line' type objects.
%
%
% 2010-04 phli
%

opts = inputParser;
opts.addParamValue('fig',        gcf);
opts.addParamValue('width',      []);
opts.addParamValue('height',     []);
opts.addParamValue('yscale',     1);
opts.addParamValue('ymode',      'global');
opts.addParamValue('tracesx',    {});
opts.addParamValue('plotargs',   {{}});
opts.addParamValue('colors',     []);
opts.addParamValue('color_prop', 'Color');
opts.parse(varargin{:});
opts = opts.Results;

width    = opts.width;
height   = opts.height;
tracesx  = opts.tracesx;
plotargs = opts.plotargs;

if isempty(width)
    dx = ipdm(positions(:,1));
    width = min(dx(dx > 0));
end

if isempty(height)
    dy = ipdm(positions(:,2));
    height = min(dy(dy > 0)) * opts.yscale;
end

if isempty(plotargs) || ~iscell(plotargs{1})
    plotargs = {plotargs};
end

num_colors = length(opts.colors);




% Set up data part 1: if traces is a simple matrix, convert it to cell array.
if isnumeric(traces)
    traces = mat2cell(traces, size(traces, 1), ones(size(traces, 2), 1));
end

% Set up data part 2
for i = 1:length(traces)
    % Put trace in column form if it's a vector
    if isvector(traces{i})
        traces{i} = traces{i}(:);
    end
    
    % Set up tracex if not given for this trace
    if i > length(tracesx) || isempty(tracesx{i})
        tracesx{i} = (1:size(traces{i}, 1))';
    end
    
    % Duplicate plotargs if not given for this trace
    if i > length(plotargs)
        plotargs{i} = plotargs{i-1};
    end
end




% Scale the y data uniformly across all traces?
if strcmp(opts.ymode, 'global')
    traceconcat = [traces{:}];
    globalmax = max(traceconcat(:));
    globalmin = min(traceconcat(:));
    clear traceconcat;
    
    yrange = globalmax - globalmin;
    yscale = height / yrange;
end


figure(opts.fig);
ax = gca;
old_hold = ishold;
hold on;
for i = 1:length(traces)
    % Rescale x into subplot position
    tracex = tracesx{i};
    centerx = positions(i,1);
    startx = centerx - width/2;
    xmin = min(tracex);
    xmax = max(tracex);
    xrange = xmax - xmin;
    xscale = width / xrange;
    xoffset = startx - (xmin*xscale);
    tracex = tracex .* xscale + xoffset;
    
    % Rescale y into subplot position
    trace = traces{i};
    centery = positions(i,2);    
    switch opts.ymode
        % Scale the y data uniformly across all traces
        case 'global'
            trace = trace .* yscale + centery;

        % Individually rescale each y
        case 'local'
            starty = centery - height/2;
            ymin = min(trace);
            ymax = max(trace);
            yrange = ymax - ymin;
            yscale = height / yrange;
            yoffset = starty - (ymin*yscale);
            trace = trace .* yscale + yoffset;
    end
    
    
    % Auto cycle through colors, if given
    plotarg = plotargs{i};
    if num_colors > 0
        icolor = mod(i-1, num_colors) + 1;
        plotarg{end+1} = opts.color_prop;
        plotarg{end+1} = opts.colors(icolor,:);
    end
    
    plot(tracex, trace, plotarg{:});
end




if ~old_hold
    hold off;
end

if nargout > 0
    varargout{1} = opts.fig;
end