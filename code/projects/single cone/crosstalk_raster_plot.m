function crosstalk_raster_plot(datarun, stimstruct, cells, varargin)

opts = inputParser;

% These are split out and passed to RASTER
opts.addParamValue('hist', true);
opts.addParamValue('hist_bin', 0.05);
opts.addParamValue('hist_color', 'b');
opts.addParamValue('hist_line_width', 2);
opts.addParamValue('color', 0.5*[1 1 1]);
opts.addParamValue('MarkerSize', 5);
opts.addParamValue('start', -0.1);
opts.addParamValue('stop', mean(diff(datarun.triggers)));

% Unmatched will also be passed along to RASTER
opts.KeepUnmatched = true; 

% These are for me; have to be added to the remove list for RASTEROPTS below.
opts.addParamValue('maps', []);

opts.parse(varargin{:});
unmatched = opts.Unmatched;
opts = opts.Results;
rasteropts = rmfield(opts, {'maps'});
rasteropts.plot = true;



intensity = unique(stimstruct.rgbs(:,1)); % These stimuli have the same intensity for R, G, B
intensity = setdiff(intensity, 0);
intensity = sort(intensity, 1, 'descend');

if isempty(opts.maps)
    maps = reshape(stimstruct.mapindices, [length(stimstruct.mapindices)/2 2])';
else
    maps = opts.maps;
end

numtrigg = length(datarun.triggers);
if numtrigg < length(stimstruct.maps)
    warning(['Not enough triggers for all pulses.  Will only use the first ' num2str(numtrigg) ' pulses']);
end
if numtrigg > length(stimstruct.maps)
    error('More triggers than maps; I don''t know how to handle this.');
end

% New format is to leave cells as cellstruct that may have empties
if iscell(cells)
    cells = cell2mat(select(cells, @(e)(~isempty(e))));
end

height = size(maps,1);
width  = size(maps,2);
for cell = cells
    figure
    
    for i = 1:size(maps,1),

        for j = 1:size(maps,2),
            map = maps(i,j);
            
            pax = subplot(height, width, sub2ind([width height], j, i));
            [res, ax] = rasterphli(datarun, cell, datarun.triggers(stimstruct.rgbs(1:numtrigg,1) == intensity & (stimstruct.maps(1:numtrigg) == map)'), rasteropts, unmatched);
            
            if ~isempty(res), numrast(j) = max(res(:,2)); end
            lax(j) = ax(1);
            rax(j) = ax(2);
            if (lax(j) > 0), set(lax(j), 'ytick', [], 'xtick', []); end
            if i ~= height, set(rax(j), 'xtick', 0); end
        end
        
        % Link axes, saving right axis ylims first
        lax = lax(lax > 0);
        if ~isempty(lax), linkaxes(lax); end;
        for j = 1:length(rax)
            ylims(j,:) = get(rax(j), 'ylim');
        end
        linkaxes(fliplr(rax));
        linkprop(fliplr(rax), 'YTickMode');
        
        if ~isempty(lax), set(lax(1), 'ylim', [0 max(numrast)]); end;
        set(rax(1), 'ylim', [min(ylims(:,1)) max(ylims(:,2))]);
        set(rax(1), 'YTickMode', 'auto');
        set(rax(1), 'xlim', [opts.start, opts.stop]);
    end
end