function [rax,histx,histy,res] = wack_raster_plot(datarun, stimstruct, cells, cellregions, varargin)
% For the old style wack, the one that was actually wack because we used
% white noise.  This was only run a few times...

opts = inputParser;

% These are split out and passed to RASTER
opts.addParamValue('hist', true);
opts.addParamValue('hist_bin', 0.05);
opts.addParamValue('hist_line_width', 2);
opts.addParamValue('color', 0.5*[1 1 1]);
opts.addParamValue('MarkerSize', 5);
opts.addParamValue('triggers', datarun.triggers);
opts.addParamValue('start', -0.1);
opts.addParamValue('stop', []);

% Unmatched will also be passed along to raster
opts.KeepUnmatched = true;

% These are mine
opts.addParamValue('link_ylims', true);
opts.addParamValue('subplot_width',  4);
opts.addParamValue('subplot_height', length(cells));
%opts.addParamValue('hist_colors', {'r', [0 0.5 0], 'b'});

opts.parse(varargin{:});
unmatched = opts.Unmatched;
opts = opts.Results;
if isempty(opts.stop)
    opts.stop = mean(diff(opts.triggers));
end;
rasteropts = keepfields(opts, 'hist', 'hist_bin', 'hist_line_width', 'color', 'MarkerSize', 'start', 'stop');


for i = 1:length(cells)
    cell = cells{i};
    if isempty(cell), continue; end
    regions = cellregions{i};
    intensities = stimstruct.intensities(regions,:);
    
    % Get trials for each condition
    trials{1} = sum(intensities) == 0; % Both on
    trials{2} = sum(intensities) == 2; % Both off
    trials{3} = intensities(1,:) == 0 & intensities(2,:) == 1; % 1 on,  2 off
    trials{4} = intensities(1,:) == 1 & intensities(2,:) == 0; % 1 off, 2 on
    
    for p = 1:4
        subplot(opts.subplot_height, opts.subplot_width, sub2ind([opts.subplot_width opts.subplot_height], p, i));
        [res{i,p}, ax, histx, histy{i,p}] = rasterphli(datarun, cell, opts.triggers(trials{p}), rasteropts, unmatched);
        if ~isempty(res{i,p}), numrast(i,p) = max(res{i,p}(:,2)); end
        lax(i,p) = ax(1);
        rax(i,p) = ax(2);
    end
    
    if opts.link_ylims
        for j = 1:length(rax(i,:))
            ylims(j,:) = get(rax(i,j), 'ylim');
        end
        linkprop(rax(i,:), 'ylim');
        set(rax(i,1), 'ylim', [min(ylims(:,1)) max(ylims(:,2))]);
    end
    
    % Set left axis ylims
    rowlax = lax(i,:);
    rowlax = rowlax(rowlax > 0);
    if ~isempty(rowlax), linkaxes(rowlax); end;
    if ~isempty(rowlax), set(rowlax(1), 'ylim', [0 max(numrast(i,:))]); end;
    for j = 1:length(rowlax),     set(rowlax(j), 'ytick', [], 'xtick', []); end
    
    % Link right axis ticks, set to auto
    linkprop(fliplr(rax(i,:)), 'YTickMode');
    set(rax(i,1), 'YTickMode', 'auto');
end

% Turn off xtick labels except for bottom row
for i = 1:(size(rax,1)-1)
    for j = 1:size(rax,2)
        if rax(i,j) > 0
            set(rax(i,j), 'XTickLabel', []);
        end
    end
end

% Turn off right tick labels except for last col
for i = 1:size(rax,1)
    for j = 1:(size(rax,2)-1)
        if rax(i,j) > 0
            set(rax(i,j), 'YTickLabel', []);
        end
    end
end


% for i = 1:size(rax,1)
%     cell = cells{i};
%     if isempty(cell), continue; end
%     subplot(opts.subplot_height, opts.subplot_width, sub2ind([opts.subplot_width opts.subplot_height], 1, i));
%     ylabel(cell);
% end