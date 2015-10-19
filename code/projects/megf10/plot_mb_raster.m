function [] = plot_mb_raster(MB, raster, trial_dur, cell_idx, cell_id, Title, x, y, print_out)
% generate subplot indices for DS raster plot with multiple dataruns. 
% So far assume that there are 8 directions.
% DS: DS structure with all datasets.
% raster: rasters of all datasets.
% print: 1:print as pdf and close
%        0:don't print

[idx, xx, yy] = subplot_idx(x, y);
ii = find(~cellfun(@isempty,raster{1}),1); % get the idx of 1st non-empty cell in raster
tpn = size(raster{1}{ii}, 2);
tt = MB{1}.theta{1}(1, :);
for time = 1:tpn
    FigHandle = figure;
%     set(FigHandle, 'Position', [1 1 1400 600])
    set(FigHandle, 'Position', [1 1 1620 1080])
%     set(FigHandle, 'Position', [1 1 1080 1080])
    for j = 1:length(MB)
        if ~isempty(raster{j}{cell_idx}) && time <= length(MB{j}.rho)
            h = subplot(xx, yy, idx(j, 1)) ; polar(tt, MB{j}.rho{time}(cell_idx, :));
            polar_theta_off(h)
            for i = 2:9
                subplot(xx, yy, idx(j, i)); plot_raster(squeeze(raster{j}{cell_idx}(1, time, i-1, :)), 0, trial_dur{j}(time))
                if i == 4
                    title(Title{j})
                end 
                if mod(idx(j, i), yy) == 1
                    ylabel('trial number')
                end
                if idx(j, i) > yy*(xx-1)
                    xlabel('time (s)')
                end
            end
        end
    end
     
    if print_out
        name = [num2str(cell_id) '_' num2str(time)];
        screen_size = [24 12];
        set(figure(1), 'paperpositionmode', 'auto');
        set(gcf, 'PaperUnits', 'inch');
        set(figure(1), 'PaperSize', screen_size);
        print(figure(1), '-dpdf', name)
        close
    end
end

        

