
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 2010-09-24-1
array_type = 519;
master_data_path = '/snle/lab/Experiments/Array/Analysis/2010-09-24-1/data021/data021';
slave_data_path = '/snle/lab/Experiments/Array/Analysis/2010-09-24-1/data022-data033/data022-data033';

master_cone_path = [single_cone_path, '2010-09-24-1_data006_data006-bayes-msf_15.00-15-test/'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load master data
clear temp_datarun
temp_datarun = load_data(master_data_path);
temp_datarun = load_params(temp_datarun,struct('verbose',1));  
temp_datarun = load_ei(temp_datarun, 'all', struct('array_type', array_type));
temp_datarun = load_sta(temp_datarun, 'load_sta', 'all', 'save_sta', false, 'keep_java_sta', false);
%temp_datarun = import_single_cone_data(temp_datarun, master_cone_path');
datarun{1} = temp_datarun;

% load slave data
clear temp_datarun
temp_datarun = load_data(slave_data_path);
temp_datarun = load_neurons(temp_datarun);
temp_datarun = load_ei(temp_datarun, 'all',struct('array_type', array_type));
datarun{2} = temp_datarun;

clear temp_datarun
threshold = 0.95; % correlaton threhold for mapping
cell_type = {4};

% EI mapping from RF run to control condition
master_cell_type = cell_type;
slave_cell_type = 'all';
[cell_list, fail_list] = map_ei(datarun{1}, datarun{2}, 'corr_threshold', threshold, 'master_cell_type', master_cell_type, 'slave_cell_type', slave_cell_type, 'verbose', true, 'troubleshoot', true);


%% -----------------------------
datarun{1} = load_sta(datarun{1}, 'load_sta', []);

datarun{1} = get_sta_fits_from_vision(datarun{1}, 'all');
plot_rf_summaries(datarun{1}, cell_type,'foa', 2, 'plot_fits', true, 'label', true, 'label_size', 10, 'fit_width', 0.25) 

print(3, '~/Desktop/off-midget-rfs.pdf', '-dpdf')

%% -----------------------------

cell_ID_A = 5701; %off midget
cell_ID_A = 6556; %off midget
cell_ID_A = 1276; %off midget
cell_ID_A = 1157; %off midget
cell_ID_A = 692; %off midget
cell_ID_A = 4996; %off midget
cell_ID_A = 5987

%cell_ID_A = 5911
%cell_ID_A = 1531
cell_ID_A = 6001

%cell_ID_A = 5191
cell_ID_A = 5476
%cell_ID_A = 6076
%cell_ID_A = 5386

cell_index_A = get_cell_indices(datarun{1}, cell_ID_A);
cell_ID_B = cell_list(cell_index_A);
cell_index_B = get_cell_indices(datarun{2}, cell_ID_B{1});

num_conditions = 12;
condition_duration = 120; % seconds
cycle_duration = 0.5; %seconds
triggers_per_cycle = 2;
bin_edges = 0:0.02:0.5;
print_flag = false;

for cnd = 1:num_conditions

    begin_time = (cnd-1) * condition_duration;
    end_time = cnd * condition_duration;

    %%% RASTER
    trigger_begin = ((begin_time ./ cycle_duration) * triggers_per_cycle) +1;
    trigger_end = (end_time ./ cycle_duration) .* triggers_per_cycle;

    cycle_trigger_indices = trigger_begin:triggers_per_cycle:trigger_end;
    num_cycles = length(cycle_trigger_indices);

    cycle_triggers = datarun{2}.triggers(cycle_trigger_indices);
    temp_spikes = datarun{2}.spikes{cell_index_B};
    temp_spike_indices = find(temp_spikes >= begin_time & temp_spikes < end_time);
    used_spikes = temp_spikes(temp_spike_indices);

    [raster_times, raster_marks] = spike_raster(used_spikes, cycle_triggers, 0, 'plot_raster', true, 'raster_size', 10);
    temp_path = ['~/Desktop/single-cone/raster-', num2str(cnd),'.pdf'];
    if print_flag
        title([num2str(cnd)])
        print(cnd, temp_path, '-dpdf')
    end
    raster_keeper{cnd} = raster_times;
    
    %%  -----------------------------------
    % make psth

    hist_matrix = zeros(length(raster_times), length(bin_edges));

    for cyc = 1:length(raster_times)
        temp_times = raster_times{cyc};
        temp_hist_times = histc(temp_times, bin_edges);
        hist_matrix(cyc,:) = temp_hist_times';
    end
    temp_psth = sum(hist_matrix)./num_cycles./bin_edges(2);
    temp_error = std(hist_matrix./bin_edges(2)) ./ sqrt(length(raster_times));
    figure(20+cnd)
    p = bar(bin_edges, temp_psth, 'histc');
    set(p, 'FaceColor', 'k')
    axis([0 0.5 0 40])
    temp_path = ['~/Desktop/single-cone/histogram-', num2str(cnd),'.pdf'];
    if print_flag
        title([num2str(cnd)])
        print(20+cnd, temp_path, '-dpdf')
    end
    
    condition_keeper(cnd).psth = temp_psth;
    condition_keeper(cnd).rate_error = temp_error;
    condition_keeper(cnd).raster = raster_times;
    
end 
%% -----------------------------
% make contrast response curves

for cnd = 1:num_conditions
    full_resp(cnd) = sum(condition_keeper(cnd).psth(3:7)) ./ 5;
    base_line(cnd) = sum(condition_keeper(cnd).psth(20:25)) ./ 6;
    resp_error(cnd) = mean(condition_keeper(cnd).rate_error(1:9));
    %full_resp(cnd) = max(condition_keeper(cnd).psth);
end
full_resp = full_resp - mean(base_line);

contrast_vals = [0.48 0.24 0.12 0.06];

figure(102); clf; 
errorbar(contrast_vals, full_resp(1:4), resp_error(1:4), 'go-');
hold on
errorbar(contrast_vals, full_resp(5:8), resp_error(1:4), 'bo-');
errorbar(2*contrast_vals, full_resp(9:12), resp_error(1:4), 'ko-');
%semilogx(2*contrast_vals, full_resp(1:4) + full_resp(5:8), 'ro-')

%semilogx(contrast_vals, full_resp(1:4) + full_resp(5:8), 'ro-')
axis([0 1 0 21])
print(102, '~/Desktop/contrast-response.pdf', '-dpdf')

regress(full_resp(1:4)', full_resp(5:8)')

regress(full_resp(1:4)',contrast_vals')
regress(full_resp(5:8)',contrast_vals')



% fit the data with an exp.

coef = [1 0];
% fit_coef_one = nlinfit(contrast_vals, full_resp(1:4), 'simple_exp_fit', coef)
% fit_coef_two = nlinfit(contrast_vals, full_resp(5:8), 'simple_exp_fit', coef)
% fit_coef_both = nlinfit(contrast_vals*2, full_resp(9:12), 'simple_exp_fit', coef)
% 
% semilogx([0.01:0.001:1], simple_exp_fit(fit_coef_one, [0.01:0.001:1]), 'g')
% semilogx([0.01:0.001:1], simple_exp_fit(fit_coef_two, [0.01:0.001:1]), 'b')
% semilogx([0.01:0.001:1], simple_exp_fit(fit_coef_both, [0.01:0.001:1]), 'k')
% hold off

x_data = [contrast_vals, contrast_vals, contrast_vals*2];
fit_coef_all = nlinfit(x_data, full_resp, 'simple_exp_fit', coef)
plot([0.01:0.001:1], simple_exp_fit(fit_coef_all, [0.01:0.001:1]), 'r')


x_data = [contrast_vals, contrast_vals, contrast_vals*2];
fit_coef_all = nlinfit(x_data, full_resp, 'square_fit', coef)
plot([0.01:0.001:1], square_fit(fit_coef_all, [0.01:0.001:1]), 'c')
hold off

% plot RF
plot_rf(datarun{1}, cell_ID_A, 'foa', 6, 'polarity', true, 'fit', false, 'scale', 5)
print(6, '~/Desktop/RF.pdf', '-dpdf')

%% ------------------------------

% load in cone files so you can figure out what was cones were stimulated
% and there respective weights.
cone_data_path = '/snle/lab/Experiments/Array/Analysis/2010-09-24-1/streamed/data006/data006';
cone_path = [single_cone_path '2010-09-24-1_data006_data006-bayes-msf_15.00-15-test'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load master data
cone_datarun = load_data(cone_data_path);
cone_datarun = load_params(cone_datarun,struct('verbose',1));  
cone_datarun = load_sta(cone_datarun, 'load_sta', 'all', 'save_sta', false, 'keep_java_sta', false);
cone_datarun = import_single_cone_data(cone_datarun, cone_path);


% FIGURES
% make figure of the streamed data
cone_datarun = get_sta_fits_from_vision(cone_datarun, 'all');
plot_rf_summaries(cone_datarun, {4},'foa', 3, 'plot_fits', true, 'label', true) 
title('streamed data')

% make figure of the analysis data
datarun{1} = get_sta_fits_from_vision(datarun{1}, 'all');
plot_rf_summaries(datarun{1}, {4},'foa', 4, 'plot_fits', true, 'label', true) 
title('analyzed data')

plot_rf_summaries(datarun{1}, {4},'foa', 3, 'plot_fits', true, 'label', false, 'fit_color','r','clear', false) 



temp_cell_id = 46;  % equals cell 5701
temp_cell_id = 83;  % equals cell 1531
temp_cell_id = 29;  % equals cell 4996


[cell_weights, cell_selection, extras] = select_cone_weights(cone_datarun, temp_cell_id,...
                                            'thresh', 0.1,...
                                            'radius', [0 inf], 'polarity', 1,...
                                            'contiguity', false,'scale', 3.0);   

cell_connectivity = cell_weights .* cell_selection; % keep weights continuous valued

[sorted_weights, sorted_indices] = sort(abs(cell_connectivity), 'descend');
sorted_weights(1:3)

weight_fraction = sorted_weights(2)./sorted_weights(1)                                        


%%
scale_up = 10;
                                        
% alternative data path
datarunC = load_data('/snle/lab/Experiments/Array/Analysis/2010-09-24-1/streamed/data006/data006');
path_and_name{1,2} = '2010-09-24-1_streamed_data006_data006-bayes-msf_15.00-15-foo';

% load data
datarunC = load_params(datarunC,struct('verbose',1));  
datarunC = load_sta(datarunC,'load_sta',[]);
datarunC = set_polarities(datarunC);
datarunC = import_single_cone_data(datarunC, path_and_name{1,2});    
datarunC.cones.mosaic = make_mosaic_struct(datarunC.cones.centers);
datarunC = set_polarities(datarunC);

datarunC = get_sta_fits_from_vision(datarunC, {4});

plot_rf_summaries(datarunC, {4},'foa', 3, 'plot_fits', true, 'label', true, 'label_size', 10, 'fit_width', 0.25) 

cell_ID = 46;

plot_rf(datarunC, cell_ID,'foa', 5, 'fit', false, 'polarity', true, 'scale', scale_up)


cell_index = get_cell_indices(datarunC, cell_ID);
datarunC = get_sta_summaries(datarunC, cell_ID);
rf = matrix_scaled_up(datarunC.stas.rfs{cell_index}, scale_up) * -1;
rf = repmat(rf, [1,1,3]);


width = datarunC.stimulus.field_width;
height = datarunC.stimulus.field_height;

%% make tessalation
centers = datarunC.cones.centers;
[V,C] = voronoin(centers);

%% make cone map
verbose = false;
clear cone_map
cone_map = zeros(width, height);
edge_map = matrix_scaled_up(cone_map, scale_up);
figure(101); clf; axis ij
new_counter = 0;
for cn = 1:length(C)
    clear temp_mask
    if all(C{cn} ~=1)
        temp_mask = poly2mask(V(C{cn},1),V(C{cn},2), height, width);
        if length(find(temp_mask ==1)) <= 40
            if verbose
                patch(V(C{cn},1),V(C{cn},2), 'r')
                axis([1 320 1 320])
                axis square
                drawnow
                pause(0.01)
            end
            new_counter = new_counter +1;
            cone_map = cone_map + (new_counter * temp_mask);
            temp_edge_map = edge(matrix_scaled_up(new_counter * temp_mask, scale_up));
            edge_map = edge_map + temp_edge_map;
            imagesc(edge_map)
            drawnow
        end
    end
end

new_edge_map = zeros([width*scale_up,height*scale_up]);
new_edge_map(edge_map>0) = 0;
edge_map(edge_map > 0) = 1;
offset_edge_map = edge_map;
offset_edge_map(offset_edge_map == 0) = NaN;

edge_map_image(1:width*scale_up,1:width*scale_up,1) = offset_edge_map;
edge_map_image(1:width*scale_up,1:width*scale_up,2) = new_edge_map;
edge_map_image(1:width*scale_up,1:width*scale_up,3) = new_edge_map;                                      

figure(61)
image(edge_map_image)
axis image
print(61, '~/Desktop/edge_map.pdf', '-dpdf')


outline_indices = find(edge_map == 1);


rf_voronoi = rf;
rf_voronoi = norm_image(rf_voronoi);

rf_red_channel = squeeze(rf_voronoi(:,:,1));
rf_green_channel = squeeze(rf_voronoi(:,:,2));
rf_blue_channel = squeeze(rf_voronoi(:,:,3));
rf_red_channel(outline_indices) = 1;
rf_green_channel(outline_indices) = 0;
rf_blue_channel(outline_indices) = 0;

rf_voronoi(:,:,1) = rf_red_channel;
rf_voronoi(:,:,2) = rf_green_channel;
rf_voronoi(:,:,3) = rf_blue_channel;

figure(1)
image(rf_voronoi)
axis image
rf_com = datarunC.stas.rf_coms{46};
rf_com = rf_com * scale_up;
window_size = 15;
begin_x = rf_com(1) - (window_size*scale_up);
end_x = rf_com(1) + (window_size*scale_up);
begin_y = rf_com(2) - (window_size*scale_up);
end_y = rf_com(2) + (window_size*scale_up);
axis([begin_x end_x begin_y end_y])



datarunD = datarunC;
% modify cone centers
datarunD.cones.centers = (datarunD.cones.centers .* scale_up) - repmat([scale_up/2-0.5, scale_up/2-0.5], length(datarunC.cones.types), 1);

plot_cone_mosaic(datarunD, 'fig_or_axes', 1,...
                'cone_size', 10, 'cone_colors', [1 1 1], 'clear', false, 'bg_color', [])
axis([begin_x end_x begin_y end_y])
           
print(1, '~/Desktop/test.pdf', '-dpdf', '-painters')

figure(2); clf;
%plot_rf(datarunC, 46, 'foa', 2, 'fit', false);
plot_cone_mosaic(datarunC, 'fig_or_axes', 2,...
    'cone_size', 4, 'cone_colors', [1 1 1], 'clear', true, 'bg_color', [0.5 0.5 0.5])

            print(2, '~/Desktop/cone-mosaic-b.pdf','-dpdf')

%%


% extract connectivity
[mosaic_weights, selection, extras] = select_cone_weights(datarunC, cell_ID,...
                                            'thresh', 0.1,...
                                            'radius', [0 inf], 'polarity', 1,...
                                            'contiguity', false,'scale', 3.0);   
connectivity = mosaic_weights .* selection; % keep weights continuous valued

% get indices to strongest cones
[num_cones, num_rgcs] = size(connectivity);
peak_cone_indices = zeros(num_rgcs,1);
for cc = 1:num_rgcs
    [max_val, temp_cn_index] = max(connectivity(:,cc));
    peak_cone_indices(cc) = temp_cn_index;
end

% make and plot the map of the strongest cones
peak_cone_map = zeros(width, height);
edge_map = matrix_scaled_up(peak_cone_map, scale_up);

verbose = true;
figure(102); clf; axis ij
new_counter = 0;
for cn = 1:length(C)
    if ~isempty(find(peak_cone_indices == cn,1))
        clear temp_mask
        if all(C{cn} ~=1)
            temp_mask = poly2mask(V(C{cn},1),V(C{cn},2), height, width);
            if length(find(temp_mask ==1)) <= 40
                if verbose
                    patch(V(C{cn},1),V(C{cn},2), 'r')
                    axis([1 320 1 320])
                    axis square
                    drawnow
                    pause(0.01)
                end
                new_counter = new_counter +1;
                peak_cone_map = peak_cone_map + temp_mask;
                if new_counter == 2
                    continue
                end
                
            end
        end
    end
end
peak_cone_map = peak_cone_map';
dlmwrite([cone_folder,'/map-0000.txt'], peak_cone_map, 'delimiter', '\t', 'newline', 'pc')

%%

% compute the spike count histograms for one condition

temp_raster = raster_keeper{1};
resp_spike_count = zeros(length(temp_raster),1);
baseline_spike_count = resp_spike_count;
for trl = 1:length(temp_raster)
    temp_times = temp_raster{trl};
    resp_indices = find(temp_times <= 0.140);
    baseline_indices = find(temp_times > 0.360);
    
    resp_spike_count(trl) = length(resp_indices);
    baseline_spike_count(trl) = length(baseline_indices);
end
temp_bins = [0,1,2,3,4,5,6];
resp_hist = hist(resp_spike_count, temp_bins);
baseline_hist = hist(baseline_spike_count, temp_bins);

plot(temp_bins, resp_hist, 'b', temp_bins, baseline_hist, 'r')

bar(temp_bins',[resp_hist;baseline_hist]')
hold on
bar(temp_bins, baseline_hist, 'k')
hold off

discriminant = mean(resp_spike_count) - mean(baseline_spike_count);
correct_indices = find(resp_spike_count >= discriminant);
A_correct = length(correct_indices) ./ length(resp_spike_count);
b_correct_indices = find(baseline_spike_count < discriminant);
B_correct = length(b_correct_indices) ./ length(baseline_spike_count);
p_correct = (A_correct + B_correct) ./2














