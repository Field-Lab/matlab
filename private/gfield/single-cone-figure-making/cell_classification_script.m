% classification

% get RF sizes


% get indices for all the cell types
cell_types = {1, 2, 3, 4, 5};
all_cell_indices = get_cell_indices(datarun, cell_types);
on_par_indices = get_cell_indices(datarun, {1});
off_par_indices = get_cell_indices(datarun, {2});
on_mid_indices = get_cell_indices(datarun, {3});
off_mid_indices = get_cell_indices(datarun, {4});
sbc_indices = get_cell_indices(datarun, {5});


num_cells = length(all_cell_indices);

temp_rf_fits = datarun.vision.sta_fits(all_cell_indices);


rf_sizes = zeros(num_cells,1);
for RGC = 1:num_cells
    temp_rf_size = geomean(temp_rf_fits{RGC}.sd)./2;
    rf_sizes(RGC) = temp_rf_size;
end
pixel_size = 5.5;
stixel_size = 5.0;
rf_sizes = rf_sizes * pixel_size * pixel_size * stixel_size;

% get TC information
summed_red_tc = zeros(num_cells, 1);
summed_green_tc = summed_red_tc;
summed_blue_tc = summed_red_tc;
for RGC = 1:num_cells

    temp_tc = datarun.stas.time_courses{all_cell_indices(RGC)};

    % get max and min of TC
    temp_max_blue = max(temp_tc(:,3));
    temp_min_blue = min(temp_tc(:,3));
    temp_max_green = max(temp_tc(:,2));
    temp_min_green = min(temp_tc(:,2));

    % get peak
    if abs(temp_max_blue) >= abs(temp_min_blue)
        temp_peak_blue = temp_max_blue;
    else
        temp_peak_blue = temp_min_blue;
    end

    if abs(temp_max_green) >=  abs(temp_min_green)
        temp_peak_green = temp_max_green;
    else
        temp_peak_green = temp_min_green;
    end

    if abs(temp_peak_green) >= abs(temp_peak_blue)
        peak_sensitivity = temp_peak_green;
    else
        peak_sensitivity = temp_peak_blue;
    end

    summed_red_tc(RGC) = sum(temp_tc(20:23,1)) ./ abs(peak_sensitivity);
    summed_green_tc(RGC) = sum(temp_tc(20:23,2)) ./ abs(peak_sensitivity);
    summed_blue_tc(RGC) = sum(temp_tc(20:23,3)) ./ abs(peak_sensitivity);

    

    if 0
        figure(1)
        clf
        hold on
        plot(temp_tc(:,1) ./ abs(peak_sensitivity), 'r')
        plot(temp_tc(:,2) ./ abs(peak_sensitivity), 'g')
        plot(temp_tc(:,3) ./ abs(peak_sensitivity), 'b')
        hold off
        pause
    end

    % PCA based classification
    red_tc(RGC,:) = temp_tc(:,1) ./ abs(peak_sensitivity);
    green_tc(RGC,:) = temp_tc(:,2) ./ abs(peak_sensitivity);
    blue_tc(RGC,:) = temp_tc(:,3) ./ abs(peak_sensitivity);

end

figure(10)
blue_signal = summed_blue_tc ./ (abs(summed_blue_tc) + abs(summed_green_tc) + abs(summed_red_tc));
green_signal = summed_green_tc ./ (abs(summed_blue_tc) + abs(summed_green_tc) + abs(summed_red_tc));
plot(rf_sizes,(green_signal), 'ko')

[junk, sbc_pointers] = intersect(all_cell_indices, sbc_indices);
[junk, on_par_pointers] = intersect(all_cell_indices, on_par_indices);
[junk, off_par_pointers] = intersect(all_cell_indices, off_par_indices);
[junk, on_mid_pointers] = intersect(all_cell_indices, on_mid_indices);
[junk, off_mid_pointers] = intersect(all_cell_indices, off_mid_indices);



figure(11)
marker_size = 12;
clf
hold on
plot(rf_sizes(sbc_pointers), blue_weights(sbc_pointers,1), 'ko', 'MarkerSize', marker_size)
plot(rf_sizes(on_par_pointers), blue_weights(on_par_pointers,1), 'ko', 'MarkerSize', marker_size)
plot(rf_sizes(off_par_pointers), blue_weights(off_par_pointers,1), 'ko', 'MarkerSize', marker_size)
plot(rf_sizes(on_mid_pointers), blue_weights(on_mid_pointers,1), 'ko', 'MarkerSize', marker_size)
plot(rf_sizes(off_mid_pointers), blue_weights(off_mid_pointers,1), 'ko', 'MarkerSize', marker_size)
xlabel('size')
ylabel('normalized blue sensitivity')
title('ganglion cell classification')
hold off    
print(11, '~/Desktop/classification','-dpdf')



% substract mean from tcs
mean_red = mean(red_tc);
mean_green = mean(green_tc);
mean_blue = mean(blue_tc);
red_tc = red_tc - repmat(mean_red, length(red_tc(:,1)), 1);
green_tc = green_tc - repmat(mean_green, length(green_tc(:,1)), 1);
blue_tc = blue_tc - repmat(mean_blue, length(blue_tc(:,1)), 1);

[red_pcs, red_weights, red_eig_vals] = princomp(red_tc);
[green_pcs, green_weights, green_eig_vals] = princomp(green_tc);
[blue_pcs, blue_weights, blue_eig_vals] = princomp(blue_tc);

figure(12)
marker_size = 12;
clf
hold on
plot(green_weights(on_par_pointers,1), rf_sizes(on_par_pointers),'k^', 'MarkerSize', marker_size)
plot(green_weights(off_par_pointers,1), rf_sizes(off_par_pointers),'kv', 'MarkerSize', marker_size)
plot(green_weights(on_mid_pointers,1), rf_sizes(on_mid_pointers),'k+', 'MarkerSize', marker_size)
plot(green_weights(off_mid_pointers,1), rf_sizes(off_mid_pointers),'kx', 'MarkerSize', marker_size)
plot(green_weights(sbc_pointers,1), rf_sizes(sbc_pointers),'ko', 'MarkerSize', marker_size)
ylabel('size')
xlabel('first PC of green time course')
title('ganglion cell classification')
axis([-2 2 0 400])
hold off
print(12, '~/Desktop/classification','-dpdf')









