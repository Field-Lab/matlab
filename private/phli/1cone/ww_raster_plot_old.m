function [rax,histx,histy,res] = ww_raster_plot_old(datarun, stimstruct, cells, varargin)

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
opts.addParamValue('hist_colors', {'r', [0 0.5 0], 'b', [1 0 1]});
opts.addParamValue('subplot_add_width', 0);
opts.addParamValue('subplot_add_height', 0);

opts.parse(varargin{:});
unmatched = opts.Unmatched;
opts = opts.Results;
if isempty(opts.stop)
    opts.stop = mean(diff(opts.triggers));
end;
rasteropts = keepfields(opts, 'hist', 'hist_bin', 'hist_line_width', 'color', 'MarkerSize', 'start', 'stop');


% Get only trials that were actually run (don't get them past number of triggers)
pulses = stimstruct.pulses(1:length(opts.triggers));


% Identify trials for different conditions
intensity = 0.48;
for i = 1:stimstruct.numcones
    up(:,i)    = cell2mat(collect(pulses, @(s) (s.rgbs{i} == num2cell(intensity * [ 1  1  1]))));
    down(:,i)  = cell2mat(collect(pulses, @(s) (s.rgbs{i} == num2cell(intensity * [-1 -1 -1]))));
    blank(:,i) = cell2mat(collect(pulses, @(s) (s.rgbs{i} == num2cell([0 0 0]))));
end
singles = sum(blank,2) == (stimstruct.numcones - 1);
doubles = sum(blank,2) == (stimstruct.numcones - 2);


subplot_width = max(stimstruct.numcones, nchoosek(stimstruct.numcones, 2)) + opts.subplot_add_width;
subplot_height = 2 + 1 + stimstruct.numcones + opts.subplot_add_height;
for i = 1:length(cells)
    cell = cells{i};
    if isempty(cell), continue; end
    
    figure;

    % Single up
    for p = 1:stimstruct.numcones
        subplot(subplot_height, subplot_width, sub2ind([subplot_width subplot_height], p, 1));
        trials = up(:,p) & singles;
        [res{1,p}, ax, histx, histy{1,p}] = rasterphli(datarun, cell, opts.triggers(trials), rasteropts, unmatched, 'hist_color', opts.hist_colors{p});
        if ~isempty(res{1,p}), numrast(1,p) = max(res{1,p}(:,2)); end
        lax(1,p) = ax(1);
        rax(1,p) = ax(2);

        title(num2str(p));
        if p == 1, ylabel('up'); end
    end
    
    % Single down
    for p = 1:stimstruct.numcones
        subplot(subplot_height, subplot_width, sub2ind([subplot_width subplot_height], p, 2));
        trials = down(:,p) & singles;
        [res{2,p}, ax, histx, histy{2,p}] = rasterphli(datarun, cell, opts.triggers(trials), rasteropts, unmatched, 'hist_color', opts.hist_colors{p});
        if ~isempty(res{2,p}), numrast(2,p) = max(res{2,p}(:,2)); end
        lax(2,p) = ax(1);
        rax(2,p) = ax(2);
        
        if p == 1, ylabel('down'); end
    end
    clear trials;
    

    % Doubles
    col = 1;
    for j = 1:stimstruct.numcones
        for k = (j+1):stimstruct.numcones
            trials{1} = up(:,j) & up(:,k) & doubles;
            trials{2} = down(:,j) & down(:,k) & doubles;
            trials{3} = up(:,j) & down(:,k) & doubles;
            trials{4} = down(:,j) & up(:,k) & doubles;
            
            for p = 1:4
                subplot(subplot_height, subplot_width, sub2ind([subplot_width subplot_height], col, p+3));
                [res{col,p}, ax, histx, histy{col,p}] = rasterphli(datarun, cell, opts.triggers(trials{p}), rasteropts, unmatched);
                if ~isempty(res{col,p}), numrast(col,p) = max(res{col,p}(:,2)); end
                lax(col,p) = ax(1);
                rax(col,p) = ax(2);
                
                if p == 1, title([num2str(j) '+' num2str(k)]); end
                if col == 1 && p == 1, ylabel('->  ->'); end
                if col == 1 && p == 2, ylabel('<-  <-'); end
                if col == 1 && p == 3, ylabel('<-  ->'); end
                if col == 1 && p == 4, ylabel('->  <-'); end
            end
            
            col = col + 1;
        end    
    end
    
end