function [maximum, minimum, minimum_ind,  maximum_ind, zc, bi_ind, lobe1, lobe2, ratio_to_neighbors]= ei_properties(datarun, cell, array_size)

% cell_list = [80 108 348 392 468 498 589];
% for i = 1:length(cell_list)
%     cell= cell_list(i);
if array_size == 512
datarun = load_ei(datarun, cell, 'array_type', 512);
    [xc,yc] = getElectrodeCoords512(); 
    cell_id = get_cell_indices(datarun, cell);

waveforms = datarun.ei.eis{cell_id};
[~,ind] = max(waveforms(:));
[max_electrode,~] = ind2sub(size(waveforms),ind);

% positions = datarun.ei.position;


difference = repmat([xc(max_electrode), yc(max_electrode)], size(yc,2),1) - [xc', yc'];
distance = [[1:size(difference,1)]', sqrt(difference(:,1).^2 + difference(:,2).^2)];
sorted_distance= sortrows(distance, 2);
index = sorted_distance(:,2) > 0 & sorted_distance(:,2) < 70;
neighbor_ids = sorted_distance(index,1);

else
    datarun = load_ei(datarun, cell, 'array_type', 519);
    [xc,yc] = getElectrodeCoords519(); 
    xc = xc';
    yc = yc';
cell_id = get_cell_indices(datarun, cell);

waveforms = datarun.ei.eis{cell_id};
[~,ind] = max(abs(waveforms(:)));
[max_electrode,~] = ind2sub(size(waveforms),ind);

% positions = datarun.ei.position;


difference = repmat([xc(max_electrode), yc(max_electrode)], size(yc,2),1) - [xc', yc'];
distance = [[1:size(difference,1)]', sqrt(difference(:,1).^2 + difference(:,2).^2)];
sorted_distance= sortrows(distance, 2);
index = sorted_distance(:,2) > 0 & sorted_distance(:,2) < 45;
neighbor_ids = sorted_distance(index,1);
end


sel_wfs = waveforms(neighbor_ids, :);
avg_height_neigh = mean(max(abs(sel_wfs')));

figure; plot(sel_wfs')
hold on
plot(waveforms(max_electrode,:), 'linewidth', 2)


max_waveform = waveforms(max_electrode, :);

    [max_ , max_ind] = max(max_waveform);
    [min_ , min_ind] = min(max_waveform);
    if max_ind < min_ind
        maximum = min_;
        minimum = max_;
        minimum_ind = max_ind;
        maximum_ind = min_ind;
    else
        maximum = max_;
        minimum = min_;
        maximum_ind = max_ind;
        minimum_ind = min_ind;
    end
    bi_ind = abs(maximum)/abs(minimum);
    [~,zc] = crossing(max_waveform);
    index = zc > minimum_ind & zc < maximum_ind;
    if sum(index) > 0 && sum(index) < 2
        try
        zc= zc(index);
        lobe1 = fwhm(1:round(zc)+1, [max_waveform(1:round(zc)), max_waveform(round(zc))*2]);
        lobe2 = fwhm(round(zc)-1:size(max_waveform,2), [max_waveform(round(zc))*2, max_waveform(round(zc):size(max_waveform,2))]);
        catch
            zc = nan;
            lobe1 = nan;
            lobe2 = nan;
        end
        else
        zc = nan;
        lobe1 = nan;
        lobe2 = nan;
    end
    ratio_to_neighbors = abs(maximum)/abs(avg_height_neigh);
% end