master = '/snle/lab/Experiments/Array/Analysis/2008-08-27-5/data005/data005/data005';
slave = '/snle/lab/Experiments/Array/Analysis/2008-08-27-5/data003/data003/data003';

master = '/snle/lab/Experiments/Array/Analysis/2008-08-27-5/data004/data004';


RGBtoLMS = [8.5667256e-7 3.2690937e-7 3.1051197e-8; 2.8230554e-6 3.051172e-6 1.5783603e-7; 6.9619523e-7 1.1124398e-6 1.8809839e-6];
RGBtoLMS = RGBtoLMS * 1e8;
LMStoRGB = inv(RGBtoLMS);
transform = RGBtoLMS';
transform2 = LMStoRGB' * 1000;
 
ConvertColor = 1;
Rebin = 1;  % whether or not to rebin datarun B
threshold = 0.3;
CellType = {4};

% load data
clear temp_datarun
temp_datarun = load_data(master);
temp_datarun = load_params(temp_datarun,struct('verbose',1));  
temp_datarun = load_ei(temp_datarun,CellType, struct('array_type', 519));
temp_datarun = load_sta(temp_datarun, 'load_sta', CellType);
datarun{1} = temp_datarun;

clear temp_datarun
temp_datarun = load_data(slave);
temp_datarun = load_params(temp_datarun,struct('verbose',1));  
temp_datarun = load_ei(temp_datarun,CellType, struct('array_type', 519));
temp_datarun = load_sta(temp_datarun, 'load_sta', CellType);
datarun{2} = temp_datarun;


[datarunA, datarunB] = map_gdf(datarun{1}, datarun{2}, 'corr_threshold', 0.95, 'master_cell_type', CellType, 'slave_cell_type', CellType, 'verbose', true, 'troubleshoot', true);

cell_indices_A = get_cell_indices(datarunA, CellType);
cell_indices_B = get_cell_indices(datarunB, CellType);

if length(cell_indices_B) ~= length(unique(cell_indices_B))
    warning('the mapping was not 1 to 1, try increasing the corr_threshold')
    length(cell_indices_B)
    length(unique(cell_indices_B))
end

% rebin the datarun B if needed
%datarunB.stas.rebinned_stas = cell(size(datarunB.stas.stas));
if Rebin
   for cll = 1:length(cell_indices_B)
       TempBinnedSTA = matrix_rebinned(datarunB.stas.stas{cell_indices_B(cll)}, 10);
       datarunB.stas.stas{cell_indices_B(cll)} = TempBinnedSTA;
   end
end

Rebin_time = 1;
time_bin_factor = 6;
if Rebin_time == 1
    for cll = 1:length(cell_indices_A)
        temp_sta = datarunA.stas.stas{cell_indices_A(cll)};
        [x_dim, y_dim, color_dim, time_dim] = size(temp_sta);
        
        num_time_bins = floor(time_dim ./ time_bin_factor);
        time_course_rem = rem(time_dim, time_bin_factor);        

        for bin = 1:num_time_bins
            % define indexing
            start_index = time_course_rem + 1 + ((bin-1) * time_bin_factor);
            end_index = time_course_rem + (bin * time_bin_factor);
            
            temp_block = sum(reshape(temp_sta(:,:,:,start_index:end_index), [], time_bin_factor), 2);
            binned_block = reshape(temp_block,x_dim,y_dim,color_dim);
            returned_sta(:,:,:,bin) = binned_block;
        end
        datarunA.stas.stas{cell_indices_A(cll)} = returned_sta;
    end
end            
            
            
% get significant stixels
stixel_params = struct('select', 'max', 'radius', 16, 'thresh', 3);
datarunA = get_sta_summaries(datarunA, CellType,'marks_params', stixel_params);

% identify the peak frame from first sta in list
set(0, 'DefaultAxesColorOrder', [1 0 0; 0 1 0; 0 0 1; 0 0 0]);
figure(1)
plot(datarunA.stas.time_courses{cell_indices_A(1)})
TempTC = abs(datarunA.stas.time_courses{cell_indices_A(1)});
[maxval, maxindex] = max(TempTC);
APeakFrame = maxindex(1)

% get the sta spatial summary for cell from A and plot
figure(2)
temp_sta = get_sta(datarunA, datarunA.cell_types{CellType{1}}.cell_ids(1));
image(norm_image(temp_sta(:,:,:,APeakFrame)))

% get the STA spatial summary for cell from B and plot
figure(3)
datarunB = get_sta_summaries(datarunB, CellType);
plot(datarunB.stas.time_courses{cell_indices_B(1)})
TempTC = abs(datarunB.stas.time_courses{cell_indices_B(1)});
[maxval, maxindex] = max(TempTC);
BPeakFrame = maxindex(1)

figure(4)
temp_sta = datarunB.stas.stas{cell_indices_B(1)};
image(norm_image(temp_sta(:,:,:,BPeakFrame)))

% put peak frame into a special holding location
datarunA.stas.peak_stixels = datarunA.stas.marks;
stixel_expand = 1;
for cll = 1:length(cell_indices_A)
    [row_num, col_num] = size(datarunA.stas.marks{cell_indices_A(cll)});
    [temp_row, temp_col] = ind2sub([row_num, col_num], (find(datarunA.stas.marks{cell_indices_A(cll)})));
    if length(temp_row) ~= 1
        length(temp_row)
        cll
        error('more than one sig stix')
    else
    new_row = (temp_row - stixel_expand):1:(temp_row + stixel_expand);
    new_col = (temp_col - stixel_expand):1:(temp_col + stixel_expand);
    new_row = repmat(new_row, 1, (2*stixel_expand));
    new_col = repmat(new_col, 1, (2*stixel_expand));
    new_row = sort(new_row);
    
    temp_roi = zeros([row_num, col_num]);
    temp_roi(new_row, new_col) = 1;
    datarunA.stas.marks{cell_indices_A(cll)} = logical(temp_roi);
    end
end


% color transform STAs
if ConvertColor
    datarunA = get_stas_color_transformed(datarunA, CellType, transform);
end

temp_datarun = datarunA;
temp_datarun.stas.stas = temp_datarun.stas.color_transformed_stas;
temp_datarun = get_sta_summaries(temp_datarun, {CellType});
datarunA.stas.color_transformed_rfs = temp_datarun.stas.rfs;

% sanity check that sta was correctly color transformed
temp_sta = datarunA.stas.color_transformed_stas{cell_indices_A(1)};
figure(5)
image(norm_image(temp_sta(:,:,:,APeakFrame)))


% Test the color transform from RGB to LMS
%datarunB = get_stas_color_transformed(datarunB, {CellType}, transform2);
% sanity check that sta was correctly color transformed
%temp_sta = datarunB.stas.color_transformed_stas{cell_indices_B(1)};
%figure(5)
%image(norm_image(temp_sta(:,:,:,BPeakFrame)))




% Regress stas
normalize_to_peak = 0;
roi_size = ((stixel_expand*2) +1)^2;
red_comp = zeros((roi_size * length(cell_indices_A)), 2);
green_comp = red_comp;
blue_comp = red_comp;
for cll = 1:length(cell_indices_A)
    temp_sta_A = datarunA.stas.color_transformed_stas{cell_indices_A(cll)};
    temp_sta_B = datarunB.stas.stas{cell_indices_B(cll)};
    
    temp_rf_A = squeeze(temp_sta_A(:,:,:,APeakFrame));
    temp_rf_B = squeeze(temp_sta_B(:,:,:,BPeakFrame));

    index_to_sig_stix = find(datarunA.stas.marks{cell_indices_A(cll)});
    
    vectorized_rf_A = [];
    vectorized_rf_B = [];
    for color_cnt = 1:length(temp_rf_A(1,1,:))
        temp_mat_A = squeeze(temp_rf_A(:,:,color_cnt));
        temp_mat_B = squeeze(temp_rf_B(:,:,color_cnt));
        vectorized_rf_A = cat(1, vectorized_rf_A, temp_mat_A(index_to_sig_stix));
        vectorized_rf_B = cat(1, vectorized_rf_B, temp_mat_B(index_to_sig_stix));
    end
    
    temp_regress_factor = vectorized_rf_A' / vectorized_rf_B';
    
    vectorized_rf_A = vectorized_rf_A ./ temp_regress_factor;
    
    if normalize_to_peak
        if stixel_expand == 1
            max_val_A = vectorized_rf_A(14);
            max_val_B = vectorized_rf_B(14);
            ave_max_val = (max_val_A + max_val_B) ./ 2;
            vectorized_rf_A = vectorized_rf_A ./ ave_max_val;
            vectorized_rf_B = vectorized_rf_B ./ ave_max_val;
        else
            error('indexing will be incorrect because code assumes a 3x3 region of interest')
        end
    end
    
    temp_start = 1 + ((cll - 1) * roi_size);
    temp_end = cll * roi_size;
    % if off parasol or off midget flip sign so that vals are positive
    if CellType{1} == 2 || CellType{1} == 4 && normalize_to_peak ~= 1
        red_comp(temp_start:temp_end,1) = -1 * vectorized_rf_A(1:roi_size);
        red_comp(temp_start:temp_end,2) = -1 * vectorized_rf_B(1:roi_size);
        green_comp(temp_start:temp_end,1) = -1 * vectorized_rf_A(1+roi_size:2*roi_size);
        green_comp(temp_start:temp_end,2) = -1 * vectorized_rf_B(1+roi_size:2*roi_size);
        blue_comp(temp_start:temp_end,1) = -1 * vectorized_rf_A(1+(2*roi_size):3*roi_size);
        blue_comp(temp_start:temp_end,2) = -1 * vectorized_rf_B(1+(2*roi_size):3*roi_size);
    end
    if CellType{1} == 1 || CellType{1} == 3 || normalize_to_peak == 1
        red_comp(temp_start:temp_end,1) = vectorized_rf_A(1:roi_size);
        red_comp(temp_start:temp_end,2) = vectorized_rf_B(1:roi_size);
        green_comp(temp_start:temp_end,1) = vectorized_rf_A(1+roi_size:2*roi_size);
        green_comp(temp_start:temp_end,2) = vectorized_rf_B(1+roi_size:2*roi_size);
        blue_comp(temp_start:temp_end,1) = vectorized_rf_A(1+(2*roi_size):3*roi_size);
        blue_comp(temp_start:temp_end,2) = vectorized_rf_B(1+(2*roi_size):3*roi_size);
    end
end



figure(1)
clf
plot(red_comp(:,1),red_comp(:,2), 'r.', 'MarkerSize', 12)
hold on
plot([-0.15 0.8],[-0.15 0.8], 'k')
axis([-0.15 0.8 -0.15 0.8])
set(gca, 'FontSize', 32, 'FontName', 'Helvetica')
axis square
print(1,'/snle/home/gfield/Desktop/red','-dpdf')
hold off


figure(2)
plot(green_comp(:,1),green_comp(:,2), 'g.', 'MarkerSize', 12)
hold on
plot([-0.35 2.5],[-0.35 2.5], 'k')
axis([-0.35 2.5 -0.35 2.5])
set(gca, 'FontSize', 32, 'FontName', 'Helvetica')
axis square
print(2,'/snle/home/gfield/Desktop/green','-dpdf')
hold off

figure(3)
plot(blue_comp(:,1),blue_comp(:,2), 'b.', 'MarkerSize', 12)
hold on
plot([-0.15 1.0],[-0.15 1.0], 'k')
axis([-0.15 1.0 -0.15 1.0])
set(gca, 'FontSize', 32, 'FontName', 'Helvetica')
axis square
print(3,'/snle/home/gfield/Desktop/blue','-dpdf')
hold off

red_scale = red_comp(:,1)' / red_comp(:,2)'
green_scale = green_comp(:,1)' / green_comp(:,2)'
blue_scale = blue_comp(:,1)' / blue_comp(:,2)'

figure(4)
plot((red_comp(:,1) ./ red_scale),red_comp(:,2), 'r.', 'MarkerSize', 12)
hold on
plot([-0.15 0.8],[-0.15 0.8], 'k')
axis([-0.15 0.8 -0.15 0.8])
set(gca, 'FontSize', 32, 'FontName', 'Helvetica')
axis square
print(4,'/snle/home/gfield/Desktop/red','-dpdf')


figure(5)
plot((green_comp(:,1) ./ green_scale),green_comp(:,2), 'g.', 'MarkerSize', 12)
hold on
plot([-0.35 2.5],[-0.35 2.5], 'k')
axis([-0.35 2.5 -0.35 2.5])
set(gca, 'FontSize', 32, 'FontName', 'Helvetica')
axis square
print(5,'/snle/home/gfield/Desktop/green','-dpdf')

figure(6)
plot((blue_comp(:,1) ./ blue_scale),blue_comp(:,2), 'b.', 'MarkerSize', 12) 
hold on
plot([-0.15 1.0],[-0.15 1.0], 'k')
axis([-0.15 1.0 -0.15 1.0])
set(gca, 'FontSize', 32, 'FontName', 'Helvetica')
axis square
print(6,'/snle/home/gfield/Desktop/blue','-dpdf')


% Regress stas
linear_plot = 1;
log_plot = 0;
normalize_to_peak = 1;
for cll = 1:length(cell_indices_A)
    temp_sta_A = datarunA.stas.color_transformed_stas{cell_indices_A(cll)};
    temp_sta_B = datarunB.stas.stas{cell_indices_B(cll)};
    
    index_to_sig_stix = find(datarunA.stas.significant_stixels{cell_indices_A(cll)});
    temp_array_size = size(datarunA.stas.significant_stixels{cell_indices_A(cll)});
    [temp_rows, temp_cols] = ind2sub(temp_array_size, index_to_sig_stix);
    
    vectorized_sta_A = reshape(temp_sta_A(temp_rows, temp_cols, :, APeakFrame), 1, []);
    vectorized_sta_B = reshape(temp_sta_B(temp_rows, temp_cols, :, BPeakFrame), 1, []);
    
    temp_regress_factor = vectorized_sta_A / vectorized_sta_B;
    
    temp_sta_A = temp_sta_A ./ temp_regress_factor;

    frame_A = temp_sta_A(:,:,:,APeakFrame);
    frame_B = temp_sta_B(:,:,:,BPeakFrame);
    
    if normalize_to_peak
        [temp_row_num, temp_col_num] = size(frame_A);
        peak_A = min(min(frame_A(:,:,2)));
        peak_B = min(min(frame_B(:,:,2)));
        peak = mean([peak_A peak_B]);
        frame_A = frame_A ./ peak;
        frame_B = frame_B ./ peak;
    end
  
    if linear_plot 
        
        figure(1)
        hold on
        plot(frame_A(temp_rows,temp_cols,2), frame_B(temp_rows,temp_cols,2), 'g.')
        plot([-1 5], [-1 5], 'k')
        axis([-0.4 1.4 -0.4 1.4])
        axis square
        
        figure(2)
        hold on
        plot(frame_A(temp_rows,temp_cols,1), frame_B(temp_rows,temp_cols,1), 'r.')
        plot([-1 5], [-1 5], 'k')
        axis([-0.4 1.4 -0.4 1.4])
        axis square
        
        figure(3)
        hold on
        plot(frame_A(temp_rows,temp_cols,3), frame_B(temp_rows,temp_cols,3), 'b.')
        plot([-1 5], [-1 5], 'k')
        axis([-0.4 1.4 -0.4 1.4])
        axis square
    end
    if log_plot
        temp_frame_A = frame_A(:,:,2);
        positive_indices = find(temp_frame_A(index_to_sig_stix) > 0);
        negative_indices = find(temp_frame_A(index_to_sig_stix) < 0);
        
        [temp_positive_rows, temp_positive_cols] = ind2sub(temp_array_size, positive_indices);
        [temp_neg_rows, temp_neg_cols] = ind2sub(temp_array_size, negative_indices);       
        
        loglog(abs(frame_A(temp_positive_rows,temp_positive_cols,1)), abs(frame_B(temp_positive_rows,temp_positive_cols)), 'r.')
        hold on
        loglog(abs(frame_A(temp_positive_rows,temp_positive_cols,2)), abs(frame_B(temp_positive_rows,temp_positive_cols,2)), 'g.')
        loglog(abs(frame_A(temp_positive_rows,temp_positive_cols,3)), abs(frame_B(temp_positive_rows,temp_positive_cols,3)), 'b.')
        
        %loglog(abs(frame_A(temp_neg_rows,temp_neg_cols, 1)), abs(frame_B(temp_neg_rows,temp_neg_cols, 1)), 'ro')
        %hold on
        %loglog(abs(frame_A(temp_neg_rows,temp_neg_cols, 2)), abs(frame_B(temp_neg_rows,temp_neg_cols, 2)), 'go')
        %loglog(abs(frame_A(temp_neg_rows,temp_neg_cols, 3)), abs(frame_B(temp_neg_rows,temp_neg_cols, 3)), 'bo')
        
        plot([0.0001 10], [0.0001 10], 'k')
    end
end
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Jeff's code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Get rfs from the STAs.
temp_sta_A = datarunA.stas.color_transformed_stas{cell_indices_A(1)};
temp_sta_B = datarunB.stas.stas{cell_indices_B(1)};
%figure(1)
%imagesc(norm_image(temp_sta(:,:,:,APeakFrame)))

temp_sig_stix = datarunA.stas.significant_stixels{cell_indices_A(1)};

temp_tc_A = time_course_from_sta(temp_sta_A, temp_sig_stix);
temp_tc_B = time_course_from_sta(temp_sta_B, temp_sig_stix);

temp_rf_A = rf_from_sta(temp_sta_A,'regress',struct('r_timecourse',temp_tc_A));
temp_rf_B = rf_from_sta(temp_sta_B,'regress',struct('r_timecourse',temp_tc_B));
%figure(2)
%imagesc(norm_image(temp_rf))

temp_com_A = rf_com(temp_rf_A);
temp_com_B = rf_com(temp_rf_B);

[x_A, y_A] = rf_profile(temp_rf_A(:,:,2), temp_com, 'radius', 5);
[x_B, y_B] = rf_profile(temp_rf_B(:,:,2), temp_com, 'radius', 5);


[average_x_A, average_y_A, error_y, bin_counts] = curve_from_binning(x_A, y_A,'bin_edges', 0:5);
[average_x_B, average_y_B, error_y, bin_counts] = curve_from_binning(x_B, y_B,'bin_edges', 0:5);

figure(3)
hold on
plot(x_A, y_A, 'g.')
plot(x_B, y_B, 'go')
hold off





