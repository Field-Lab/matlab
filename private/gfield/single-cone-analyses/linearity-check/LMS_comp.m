% check that large stixel and small stixel LMS are consistent

master = '/Analysis/gfield/2008-08-27-5/data005/data005';
slave = '/Analysis/gfield/2008-08-27-5/data007/data007';

master_cell_type = {3};
slave_cell_type = 'all';

ConvertColor = false;
slave_rebin = true;  % whether or not to rebin datarun B
master_rebin = false
threshold = 0.3;

% color transform info
RGBtoLMS = [8.5667256e-7 3.2690937e-7 3.1051197e-8; 2.8230554e-6 3.051172e-6 1.5783603e-7; 6.9619523e-7 1.1124398e-6 1.8809839e-6];
RGBtoLMS = RGBtoLMS * 1e8;
LMStoRGB = inv(RGBtoLMS);
transform = RGBtoLMS';

% load master datarun
clear temp_datarun
temp_datarun = load_data(master);
temp_datarun = load_params(temp_datarun,struct('verbose',1));  
temp_datarun = load_ei(temp_datarun, master_cell_type, struct('array_type', 519));
temp_datarun = load_sta(temp_datarun, 'load_sta', master_cell_type);
datarun{1} = temp_datarun;

% load slave datarun
clear temp_datarun
temp_datarun = load_data(slave);
temp_datarun = load_params(temp_datarun,struct('verbose',1));  
temp_datarun = load_ei(temp_datarun,slave_cell_type, struct('array_type', 519));
temp_datarun = load_sta(temp_datarun, 'load_sta', slave_cell_type);
datarun{2} = temp_datarun;


% map
[datarunA, datarunB] = map_gdf(datarun{1}, datarun{2}, 'corr_threshold', 0.95, 'master_cell_type', master_cell_type, 'slave_cell_type', {'ON'}, 'verbose', true, 'troubleshoot', true);

cell_indices_A = get_cell_indices(datarunA, master_cell_type);
cell_indices_B = get_cell_indices(datarunB, master_cell_type);

if length(cell_indices_B) ~= length(unique(cell_indices_B))
    warning('the mapping was not 1 to 1, try increasing the corr_threshold')
    length(cell_indices_B)
    length(unique(cell_indices_B))
end

Rebin = true;
if Rebin
   for cll = 1:length(cell_indices_A)
       TempBinnedSTA = matrix_rebinned(datarunA.stas.stas{cell_indices_A(cll)}, 2);
       datarunA.stas.stas{cell_indices_A(cll)} = TempBinnedSTA;
   end
end

% get significant stixels
stixel_params = struct('select', 'max', 'radius', 16, 'thresh', 3);
datarunA = get_sta_summaries(datarunA, master_cell_type,'marks_params', stixel_params);

% identify the peak frame from first sta in list
set(0, 'DefaultAxesColorOrder', [1 0 0; 0 1 0; 0 0 1; 0 0 0]);
figure(2)
plot(datarunA.stas.time_courses{cell_indices_A(1)})
TempTC = abs(datarunA.stas.time_courses{cell_indices_A(1)});
[maxval, maxindex] = max(TempTC);
APeakFrame = maxindex(1)

% get the sta spatial summary for cell from A and plot
figure(1)
temp_sta = get_sta(datarunA, datarunA.cell_types{3}.cell_ids(1));
image(norm_image(temp_sta(:,:,:,APeakFrame)))

% get the STA spatial summary for cell from B and plot
datarunB = get_sta_summaries(datarunB, master_cell_type);
%BTimeCourse = time_course_from_sta(datarunB.stas.stas{cell_indices_B(1)}, datarunA.stas.significant_stixels{cell_indices_A(1)});
figure(3)
plot(datarunB.stas.time_courses{cell_indices_B(1)})
TempTC = abs(datarunB.stas.time_courses{cell_indices_B(1)});
[maxval, maxindex] = max(TempTC);
BPeakFrame = maxindex(1)

figure(4)
temp_sta = datarunB.stas.stas{cell_indices_B(1)};
image(norm_image(temp_sta(:,:,:,BPeakFrame)))

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

% Regress stas
normalize_to_peak = 1;
roi_size = ((stixel_expand*2) +1)^2;
red_comp = zeros((roi_size * length(cell_indices_A)), 2);
green_comp = red_comp;
blue_comp = red_comp;
for cll = 1:length(cell_indices_A)
    temp_sta_A = datarunA.stas.stas{cell_indices_A(cll)};
    temp_sta_B = datarunB.stas.stas{cell_indices_B(cll)};
    
    temp_rf_A = squeeze(temp_sta_A(:,:,:,APeakFrame));
    temp_rf_B = squeeze(temp_sta_B(:,:,:,BPeakFrame));

    index_to_sig_stix = find(datarunA.stas.marks{cell_indices_A(cll)});
    %temp_array_size = size(datarunA.stas.marks{cell_indices_A(cll)});
    %[temp_rows, temp_cols] = ind2sub(temp_array_size, index_to_sig_stix);    
    
    vectorized_rf_A = [];
    vectorized_rf_B = [];
    for color_cnt = 1:length(temp_rf_A(1,1,:))
        temp_mat_A = squeeze(temp_rf_A(:,:,color_cnt));
        temp_mat_B = squeeze(temp_rf_B(:,:,color_cnt));
        vectorized_rf_A = cat(1, vectorized_rf_A, temp_mat_A(index_to_sig_stix));
        vectorized_rf_B = cat(1, vectorized_rf_B, temp_mat_B(index_to_sig_stix));
    end
    
    %vectorized_sta_A = reshape(temp_sta_A(temp_rows, temp_cols, :, APeakFrame), 1, []);
    %vectorized_sta_B = reshape(temp_sta_B(temp_rows, temp_cols, :, BPeakFrame), 1, []);
    
    temp_regress_factor = vectorized_rf_A' / vectorized_rf_B';
    
    %temp_sta_A = temp_sta_A ./ temp_regress_factor;
    vectorized_rf_A = vectorized_rf_A ./ temp_regress_factor;
    
    %frame_A = squeeze(temp_sta_A(:,:,:,APeakFrame));
    %frame_B = squeeze(temp_sta_B(:,:,:,BPeakFrame));
    
    %if normalize_to_peak
    %    [temp_row_num, temp_col_num] = size(frame_A);
    %    peak_A = min(min(frame_A(:,:,2)));
    %    peak_B = min(min(frame_B(:,:,2)));
    %    peak = mean([peak_A peak_B]);
    %    frame_A = frame_A ./ peak;
    %    frame_B = frame_B ./ peak;
    %end

    temp_start = 1 + ((cll - 1) * roi_size);
    temp_end = cll * roi_size;
    red_comp(temp_start:temp_end,1)  = vectorized_rf_A(1:roi_size);
    red_comp(temp_start:temp_end,2)  = vectorized_rf_B(1:roi_size);
    green_comp(temp_start:temp_end,1)  = vectorized_rf_A(1+roi_size:2*roi_size);
    green_comp(temp_start:temp_end,2)  = vectorized_rf_B(1+roi_size:2*roi_size);
    blue_comp(temp_start:temp_end,1)  = vectorized_rf_A(1+(2*roi_size):3*roi_size);
    blue_comp(temp_start:temp_end,2)  = vectorized_rf_B(1+(2*roi_size):3*roi_size);
    
end

figure(1)
plot(red_comp(:,1),red_comp(:,2), 'r.', [-0.03 0.12],[-0.03 0.12], 'k')
figure(2)
plot(green_comp(:,1),green_comp(:,2), 'g.', [-0.03 0.12],[-0.03 0.12], 'k')


    
    