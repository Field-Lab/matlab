%% define various quantitites
cone_weight_params.thresh = 0.03;
cone_weight_params.radius = [0 4];
cone_weight_params.polarity = 0;
cone_weight_params.contiguity = false;
cone_weight_params.scale = 3.0;

make_figure_flag = true;
matrix_scale_factor = 10;

dataname = 'plantain';

%% load data
[datarun_junk, cone_info] = load_data_and_cones(dataname, 'sta_summaries', false);

%% load datarun from server
cd([single_cone_path, 'saved/']);
load(dataname)

%% Get the STC kernals for the cell of interest

rgc_tcs = cone_info.rgc_tcs;
spike_times = cone_info.spike_times;
cone_inputs = cone_info.cone_inputs;
cone_ids = cone_info.cone_ids;
rgc_ids = cone_info.rgc_ids;
spike_rate = cone_info.spike_rate;
Wc = cone_info.Wc;
height = datarun.stimulus.field_height;
width = datarun.stimulus.field_width;
datarun.cones = datarun_junk.cones;
datarun = set_polarities(datarun);
clear cone_info datarun_junk


tc_length = size(rgc_tcs,2);

%cell_id = 286; % good off parasol from plantain
%cell_id = 4351; % off parasol from plantain

%cell_id = 1426; % good off midget from plantain
%cell_id = 376; % plantain off midget
%cell_id = 766; % plantain off midget
cell_id = 976; % Amazing (!!!) plantain off midget
%cell_id = 766; % off midget


%cell_id = 1486; %peach, midget 
%cell_id = 3334; % peach parasol
%cell_id = 3032 % peach, off midget

%cell_id = 302; % apple midget
%cell_id = 1369; %apple, on parasol 
%cell_id = 1081; % apple on parasol !!!

cell_index = get_cell_indices(datarun, cell_id);

window_size = 20;
temp_com = round(datarun.stas.rf_coms{cell_index,:});
bg_window_x = (temp_com(1) - window_size); %* plot_matrix_scale_factor;
ed_window_x = (temp_com(1) + window_size); %* plot_matrix_scale_factor;
bg_window_y = (temp_com(2) - window_size); %* plot_matrix_scale_factor;
ed_window_y = (temp_com(2) + window_size); %* plot_matrix_scale_factor;

% COMPUTE THE SNLS IF NEEDED
if 0
    start_stim = floor(1+start_time/(datarun.stimulus.java_movie.getRefreshTime/1000));
    end_stim = floor(1+end_time/(datarun.stimulus.java_movie.getRefreshTime/1000));
    refresh_time = datarun.stimulus.java_movie.getRefreshTime/1000;
    datarun = get_snls(datarun, cell_id, 'start_stim', start_stim, 'end_stim', end_stim./10);
    snl_struct = datarun.stas.snls{cell_index};
    plot_snl_(snl_struct.gen_signal, snl_struct.spikes, 'fit', snl_struct.fit_params, 'foa', 1)
end

cell_spike_times = spike_times{cell_index};



%%
% --- GET STA AND STC IN CONE SPACE ---

% determine significant cones for the given rgc
[mosaic_weights, selection, extras] = select_cone_weights(datarun, cell_id, cone_weight_params);
temp_cone_indices = find(selection);


% get time course from cell of interest, zero the mean and make unit lenth
temp_tc = rgc_tcs(cell_index,:);
temp_tc = temp_tc - mean(temp_tc);
temp_tc = temp_tc ./ norm(temp_tc);

% flip the impulse response to estimate the temporal filter
reverse_indices = tc_length:-1:1;
impulse_filter = temp_tc(reverse_indices);

% filter the cone projected stimulus by the RGC time course
filtered_cone_inputs = filter(impulse_filter, 1, cone_inputs(:,temp_cone_indices));



% --- COMSTRUCT THE STE AND COMPUTE STA AND STC ---


% get the covariance matrix of the stimulus (will be subtracted from the STC)
stimulus_cov = cov(filtered_cone_inputs);

% CONSTRUCT THE SPIKE-TRIGGERED ENSEMBLE
STE = filtered_cone_inputs(spike_times{cell_index},:);


% --- COMPUTE STA ---
temp_polarity = datarun.stas.polarities{cell_index};
STA = mean(STE);
STA = STA ./ max(STA) ./2;
STA = STA * temp_polarity;


% --- PLOT STA ---
figure(1)
cone_STA = Wc(:,temp_cone_indices) * STA';
cone_rf = reshape(cone_STA, [height, width, 3]);
cone_rf = cone_rf(bg_window_y:ed_window_y,bg_window_x:ed_window_x,:);
if make_figure_flag
    cone_rf = matrix_scaled_up(cone_rf, matrix_scale_factor);
end
imagesc(norm_image(cone_rf));
axis off; axis image
title('cone space from STA')
if make_figure_flag
    print(1, '~/Desktop/STA.pdf', '-dpdf')
end

%imshow(rgb2gray(norm_image(cone_rf)))

% compute and plot the raw RF for the cell
datarun = get_sta_summaries(datarun, cell_id);
temp_rf = get_rf(datarun, cell_id);
temp_rf = temp_rf(bg_window_y: ed_window_y, bg_window_x:ed_window_x,:);
temp_rf = temp_rf * datarun.stas.polarities{get_cell_indices(datarun, cell_id)};
if make_figure_flag
    temp_rf = matrix_scaled_up(temp_rf, matrix_scale_factor);
end
figure(51)
imagesc(norm_image(temp_rf))
axis off; axis square
if make_figure_flag
    print(51, '~/Desktop/raw_sta.pdf', '-dpdf')
end



% --- COMPUTE STC ---
STA_ = mean(STE)';
z_STE = STE';
z_STE = z_STE - STA_*(STA_'*z_STE)./ norm(STA_)^2;  % project out the mean
z_filtered_cone_inputs = filtered_cone_inputs';
z_filtered_cone_inputs = z_filtered_cone_inputs - STA_*(STA_'*z_filtered_cone_inputs)./ norm(STA_)^2;  % project out the mean

STC = cov(z_STE') - cov(z_filtered_cone_inputs');

% --- FACTORIZE THE STC ---
[PCs, eig_vals_matrix] = eig(STC);
eig_vals = diag(eig_vals_matrix);
[eig_vals, sorted_eig_indices] = sort(eig_vals, 'descend');
PCs = PCs(:, sorted_eig_indices);

% get indices to meaningful dimensions
keep_indices = find(abs(eig_vals) > 1e-10); % note this number might need to be changed depending on data
% store eigvals;
eig_vals = eig_vals(keep_indices); % ./ abs(sum(s_eig_vals));
PCs = PCs(:,keep_indices);

% --- PLOT THE SORTED EIGENVALUE SPECTRUM ---
figure(2)
plot(eig_vals, 'ko')

% --- PLOT EIGENVECTORS OF THE STC with increased variance---
num_dims = 10;
for dm = 1:num_dims
    image_PC = Wc(:,temp_cone_indices) * PCs(:,dm);
    image_PC = reshape(image_PC, [height, width, 3]);
    image_PC = image_PC(bg_window_y:ed_window_y,bg_window_x:ed_window_x,:);
    figure(dm+2)    
    if make_figure_flag
        image_PC = matrix_scaled_up(image_PC, matrix_scale_factor);
    end
    image(norm_image(image_PC))
    axis off; axis square
    title(['D',num2str(dm)])
    drawnow
    print_path = ['~/Desktop/',num2str(dm+2),'.pdf'];
    if make_figure_flag
        print(dm+2, print_path, '-dpdf' ,'-painters')
    end
end

% make a histogram of the mean cone weights across STC kernels
if 1
    figure(52)
    mean_cone_weights = mean(PCs);
    std_mean_cone_weights = std(mean_cone_weights);
    mean_cone_weights = mean_cone_weights ./ std_mean_cone_weights;
    hist(mean_cone_weights, [-6:0.5:5]);
    axis([-6 6 0 25])
    if make_figure_flag
        print(52, '~/Desktop/weight_histogram.pdf','-dpdf')
    end
end

%% Significance testings
num_sig_dims = 7;

[mean_spectrum, sd_spectrum] = test_stc_significance(cell_spike_times,...
                                    filtered_cone_inputs, num_sig_dims,...
                                    'iternation_number', 50);


num_sds = 4;
high_shifted_egvals = mean_spectrum + num_sds*sd_spectrum;
low_shifted_egvals = mean_spectrum - num_sds*sd_spectrum;

figure(103); clf;
max_eig_val = max(eig_vals);
plot(eig_vals./max_eig_val, 'k.', 'MarkerSize', 20)
hold on
temp_counter = (num_sig_dims+1):1:length(eig_vals);
plot(temp_counter,high_shifted_egvals./max_eig_val, 'r-', 'LineWidth', 2)
plot(temp_counter,low_shifted_egvals./max_eig_val, 'r-', 'LineWidth', 2)
axis([0, size(PCs,1)+1, -0.5 1.1])
title(['PC spectrum ', num2str(num_sig_dims)])
hold off
print(103, '~/Desktop/spectrum.pdf', '-dpdf')

%% construct kernel matrix, sta is first entry
dims_to_keep = num_sig_dims;
kernel_matrix = zeros(dims_to_keep+1, size(PCs,1));
kernel_matrix(1,:) = STA;
for dm = 1:dims_to_keep
    kernel_matrix(dm+1,:) = PCs(:,dm);
end



%% get static nonlinearity for each kernel.
verbose_plot = true;
clear kernel_params
for dm = 1:dims_to_keep+1
    
    temp_kernel_generator = kernel_matrix(dm,:) *filtered_cone_inputs' * datarun.stas.polarities{cell_index};
    cell_spike_rate = spike_rate(cell_index,:);

    if dm == 1
        [fit_params, binned_spike_rate, bin_vals] = fit_spike_response(temp_kernel_generator,...
                            cell_spike_rate, 'verbose', verbose_plot, 'fig_num', dm, 'fit_function', 'exp');
    else    
        [fit_params, binned_spike_rate, bin_vals] = fit_spike_response(temp_kernel_generator,...
                    cell_spike_rate, 'verbose', verbose_plot, 'fig_num', dm, 'fit_function', 'parabola');
    end

    kernel_params(dm).fit_params = fit_params;
    kernel_params(dm).binned_spike_rate = binned_spike_rate;
    kernel_params(dm).bin_vals = bin_vals;
    kernel_params(dm).kernel_generator = temp_kernel_generator;
    kernel_params(dm).kernel = kernel_matrix(dm,:);
    
    file_path_name = ['~/Desktop/',num2str(dm),'.pdf'];
    %print(dm, file_path_name, '-dpdf')
    
end

%% Look at joint distribution of 2 kernels


%Generate grid for 2-D projected view of intensities: 
num_bins = 10;
xb = linspace(-0.75,0.75, num_bins+1);
delta_x = xb(2) - xb(1);

joint_spike_density = zeros(num_bins);
joint_stim_density = zeros(num_bins);
a_marginal = zeros(num_bins,1);
b_marginal = zeros(num_bins,1);

kernel_a_generator = kernel_params(2).kernel_generator;
kernel_b_generator = kernel_params(3).kernel_generator;

for i_bn = 1:num_bins
    
    % get indices for generator signals falling in the ith bin
    kernel_a_indices = find(kernel_a_generator > xb(i_bn) & kernel_a_generator < xb(i_bn+1));
    a_marginal(i_bn) = length(find(cell_spike_rate(kernel_a_indices) > 0)) ./ length(kernel_a_indices);

    for j_bn = 1:num_bins
        
        % get indices for generator signals falling in the jth bin
        kernel_b_indices = find(kernel_b_generator > xb(j_bn) & kernel_b_generator < xb(j_bn + 1));
        
        if i_bn == 1
            b_marginal(j_bn) = length(find(cell_spike_rate(kernel_b_indices) > 0)) ./ length(kernel_b_indices);
        end
        
        % get indices that fall in the i-jth bin
        joint_indices = intersect(kernel_a_indices, kernel_b_indices);
        
        if ~isempty(joint_indices)
            % get spike number in this bin
            temp_spike_number = length(find(cell_spike_rate(joint_indices) > 0));

            % normalized by number of frames to get spike rate in that bin
            temp_spike_rate = temp_spike_number ./ length(joint_indices);
            
            % store value in joint_spike_densit
            joint_spike_density(i_bn, j_bn) = temp_spike_rate;
            joint_stim_density(i_bn,j_bn) = length(joint_indices);
            
        else
            
            % if there were not events in this bin, assign NaN
            joint_spike_density(i_bn, j_bn) = NaN;
        end        
    end
end

% try plotting the joint spike density in a few different ways.
figure(1)
pcolor(xb(1:num_bins),xb(1:num_bins),joint_spike_density)
colormap copper
axis square

%pcolor(xb(1:num_bins),xb(1:num_bins),joint_stim_density)

figure(2)
contour(xb(1:num_bins),xb(1:num_bins),joint_spike_density, [0.1 0.125  0.15 0.2 0.3, 0.35])
% normalize the densities


%% compute the mutual information between two kernel outputs
temp_info = 0;
for i_bn = 1:num_bins
    for j_bn = 1:num_bins
        if isnan(joint_spike_density(i_bn,j_bn))
            temp_info = temp_info +0;
        else
        temp_info = joint_spike_density(i_bn,j_bn)...
             * log2(joint_spike_density(i_bn,j_bn) ./ (a_marginal(i_bn) * b_marginal(j_bn))) * delta_x^2;
        end
    end
end
temp_info




%% compute information about stimulus from kernels
num_bins = 40;
info_spectrum = zeros(dims_to_keep,1);
for dm = 1:dims_to_keep
    kernel = kernel_params(dm).kernel;
    
    stimulus_dist = filtered_cone_inputs * kernel';
    
    frames_with_spikes = find(cell_spike_rate > 0);
    
    spike_triggered_dist = filtered_cone_inputs(frames_with_spikes,:) * kernel';
    
    [stim_hist, temp_bins] = hist(stimulus_dist, num_bins);
    [spike_trig_hist, temp_bins] = hist(spike_triggered_dist, temp_bins);
    delta_bin = temp_bins(2) - temp_bins(1);
    
%     figure(1); clf;
%     bar(temp_bins, stim_hist, 'k')
%     hold on
%     bar(temp_bins, spike_trig_hist,'r')
%     hold off

    % normalize
    stim_dist = stim_hist ./ length(stimulus_dist);
    spike_trig_dist = spike_trig_hist ./ length(spike_triggered_dist);
    
    % calculate the information
    temp_info = 0;
    for bn = 1:num_bins
        if stim_dist(bn) == 0 || spike_trig_dist(bn) == 0
            temp_info = temp_info + 0;
        else
            %spike_dist_entropy = spike_trig_dist(bn) * log2(spike_trig_dist(bn) ./  stim_dist(bn)) * delta_bin;                                   ) * delta_bin;
            %stim_dist_entropy = stim_dist(bn) * log2(stim_dist(bn)) * delta_bin;
            temp_info = temp_info + (spike_trig_dist(bn) * log2(spike_trig_dist(bn) ./  stim_dist(bn)));
    
        end
        clear
            
    end
    info_spectrum(dm)  = temp_info;
    
    
end
    
%% Compute the information in the joint distributions




%% build the model neuron

% % LNP neuron
% NL_output = exp_fit(kernel_params(1).fit_params, kernel_params(1).kernel_generator);
% NL_output(NL_output <0) = 0;
% LNP_spikes = poissrnd(NL_output);
% LNP_error = sqrt(mean((cell_spike_rate - LNP_spikes).^2));
% 
% plot(cell_spike_rate(1:1000), 'r');
% hold on
% plot(LNP_spikes(1:1000), 'b')
% hold off
% 
% 
% 
% % MDLP
% % split data in half to fit
% stop_point = floor(length(cell_spike_rate) * 0.5);
% fit_range = [1,stop_point];
% validate_range = [stop_point+1,length(cell_spike_rate)];
% 
% weights = ones(1,length(kernel_params));
% weights(3) = 0.5;
% fit_fun = @(weights)MDLP_error(weights, cell_spike_rate, kernel_params, fit_range);
% 
% [fit_weights, error_val] = fminsearch(fit_fun, weights);
% fit_weights
% error_val
% 
% 
% fit_error = MDLP_error(fit_weights, cell_spike_rate, kernel_params, validate_range)
% 
% 
% 
% 
% %%
% exp_fit_params = kernel_params(1).fit_params;
% %%
% 
% 
% 
% num_iters = 100;
% temp_tolerance = 0.001;
% for iter = 1:num_iters
%     clear regress_weights
%     generator_sigs = zeros(length(kernel_params),length(cell_spike_rate));
%     %generator_sigs = zeros(1,length(cell_spike_rate));
%     for dm = 1:length(kernel_params)
%     %for dm = 1:1
%         if dm == 1
%             temp_sigs= kernel_params(dm).kernel_generator;
% 
%             %generator_sigs(dm,:) = exp_fit(exp_fit_params, temp_sigs);
%             temp_sigs(temp_sigs<0) = 0;
%             generator_sigs(dm,:) = temp_sigs.^2;
%         else
%             generator_sigs(dm,:) = kernel_params(dm).kernel_generator.^2;
%             %generator_sigs(dm,:) = parabola_fit(kernel_params(dm).fit_params, kernel_params(dm).kernel_generator);
%         end
%     end
%     
%     %STA_mse = sqrt(mean((generator_sigs(1,:) - cell_spike_rate).^2))
% 
%     regress_weights = cell_spike_rate / generator_sigs;
% 
%     total_gen_sig = regress_weights * generator_sigs;
% 
%     full_model_mse = sqrt(mean((cell_spike_rate - total_gen_sig).^2))
% 
%     
%     % take ouput of MDNP and predice spike rate.  Get error.  Minimzie this
%     % error with respect to NL params
% 
% 
%     fit_fun = @(exp_params)MDLP_error(regress_weights, exp_params, cell_spike_rate, kernel_params, []);
% 
%     [fit_exp_params, error_val] = fminsearch(fit_fun, exp_fit_params);
% 
%     diff_params = fit_exp_params - exp_fit_params;
%     
%     storred_norm_diff(iter) = norm(diff_params);
%    
%     if norm(diff_params) < temp_tolerance
%         break
%     end
%     
%     exp_fit_params = fit_exp_params;
%     
% end
% iter
% norm(diff_params)
% 
% STA_mse
% full_model_mse
% 
% num_pts_per_bin = 500;
% [sorted_generator_vals, sorted_generator_indices] = sort(total_gen_sig,'ascend');
% num_bins = floor(length(total_gen_sig)./num_pts_per_bin);
% bin_vals = zeros(num_bins,1);
% binned_spike_rate = zeros(num_bins,1);
% for bn = 1:num_bins
%     begin_bin = 1 + ((bn-1) * num_pts_per_bin);
%     end_bin = ((bn-1) * num_pts_per_bin) + num_pts_per_bin;
%     temp_gen_indices = sorted_generator_indices(begin_bin:end_bin);
%     temp_gen_signals = sorted_generator_vals(begin_bin:end_bin);
%     
%     mean_of_bin = mean(temp_gen_signals);
%     mean_spike_rate_in_bin = mean(cell_spike_rate(temp_gen_indices));
%     
%     bin_vals(bn) = mean_of_bin;
%     binned_spike_rate(bn) = mean_spike_rate_in_bin;
% end
% 
% 
% figure
% plot(bin_vals, binned_spike_rate, 'r.-');
% hold on
% plot(bin_vals, bin_vals, 'b')
% hold off
% 
% xlabel('predicted spike rate')
% ylabel('actual spike count')
% %title('parabola fit w/ reoptimization') 
% %legend('data', 'binned data', 'fit')
% %legend boxoff


%% --------------------------------------


kernel_generators = zeros(length(kernel_params), length(cell_spike_rate));
for dm = 1:length(kernel_params)
    kernel_generators(dm,:) = kernel_params(dm).kernel_generator;
end

[num_dims, num_frames] = size(kernel_generators);

weights = ones(1, num_dims);
x_offset = 0;
%y_offset = 0;
%nl_y_offset = 0;
%offsets = [x_offset, y_offset, nl_y_offset];
offsets = x_offset;

fit_fun = @(offsets)multi_kernel_fit_error(weights, offsets, kernel_generators, cell_spike_rate);
%fit_fun = @(offsets)multi_kernel_fit_error(1, offsets, kernel_generators(1,:), cell_spike_rate);


[fit_offsets, fit_error] = fminsearch(fit_fun, offsets);

kernel_outputs = multi_kernel_output(fit_offsets, kernel_generators);
regressed_weights = cell_spike_rate / kernel_outputs


%--------
old_error = 1000
for iter = 1:100
    % this section should be run repeatedly to get a stable fit.
    fit_fun = @(offsets)multi_kernel_fit_error(regressed_weights, offsets, kernel_generators, cell_spike_rate);
    [fit_offsets, fit_error] = fminsearch(fit_fun, fit_offsets);


    kernel_outputs = multi_kernel_output(fit_offsets, kernel_generators);
    regressed_weights = cell_spike_rate / kernel_outputs;

    total_gen_sig = regressed_weights * kernel_outputs;

    temp_error = sqrt(mean((total_gen_sig - cell_spike_rate).^2));
    
    if (old_error - temp_error) < 0.0001
        break
    else
        old_error = temp_error
    end
    
end


%------




num_pts_per_bin = 500;
[sorted_generator_vals, sorted_generator_indices] = sort(total_gen_sig,'ascend');
num_bins = floor(length(total_gen_sig)./num_pts_per_bin);
bin_vals = zeros(num_bins,1);
binned_spike_rate = zeros(num_bins,1);
for bn = 1:num_bins
    begin_bin = 1 + ((bn-1) * num_pts_per_bin);
    end_bin = ((bn-1) * num_pts_per_bin) + num_pts_per_bin;
    temp_gen_indices = sorted_generator_indices(begin_bin:end_bin);
    temp_gen_signals = sorted_generator_vals(begin_bin:end_bin);
    
    mean_of_bin = mean(temp_gen_signals);
    mean_spike_rate_in_bin = mean(cell_spike_rate(temp_gen_indices));
    
    bin_vals(bn) = mean_of_bin;
    binned_spike_rate(bn) = mean_spike_rate_in_bin;
end


figure
plot(bin_vals, binned_spike_rate, 'r.-');
hold on
plot(bin_vals, bin_vals, 'b')
hold off

xlabel('predicted spike rate')
ylabel('actual spike count')
%title('parabola fit w/ reoptimization') 
%legend('data', 'binned data', 'fit')
%legend boxoff

% plot weight spectrum
figure(31)
plot([0:1:num_sig_dims], regressed_weights, 'ko-')
axis([0 num_sig_dims 0 3.0])
print(31, '~/Desktop/weights.pdf', '-dpdf')


%% With the model fit above, generate a prediction for CRGs

% make a contrast reversing stimulus
x_vals = 1:1:320; % stimulus width
y_vals = 1:1:320; % stimulus height
t_vals = 1:1:120; % stimulus duration

spatial_period = [128 96 64 48 32 24 20 16 14 12 8 4];
%spatial_period = 16;
temporal_period = 30;

mean_f_one = zeros(length(spatial_period),1);
mean_f_two = zeros(length(spatial_period),1);
sd_f_one = mean_f_one;
sd_f_two = mean_f_two;
max_f_one = mean_f_one;

x = [0:6:((length(impulse_filter)-1)*6)];
xi = [0:1:((length(impulse_filter)-1)*6)];
upsampled_inpulse_filter = interp1(x, impulse_filter,xi);


for sf = 1:length(spatial_period);
    
    spatial_phase = spatial_period(sf) .* [1,2,3,4,5,6,7,8] ./ 9 ;
    f_one = zeros(length(spatial_phase),1);
    f_two = f_one;
    
    for sp = 1:length(spatial_phase)
        
        x_amplitude = sin((2*pi*x_vals./spatial_period(sf)) + spatial_phase(sp));
        cgr_movie = repmat(x_amplitude, [length(y_vals), 1, 3, length(t_vals)]);
        
        figure(3)
        imagesc(norm_image(-1*cgr_movie(:,:,:,1)))
        
        t_amplitude = sin((2*pi*t_vals./temporal_period));
        
        for fm = 1:length(t_vals);
            cgr_movie(:,:,:,fm) = 0.5*((cgr_movie(:,:,:,fm) .* t_amplitude(fm))+1);
            %tmp_frame = 0.5*((cgr_movie(:,:,:,fm) .* t_amplitude(fm))+1);
            %test_movie(fm) = im2frame(tmp_frame);
            image(cgr_movie(:,:,:,fm))
            F(fm) = getframe;
        end
        
        %movie(F,1)
        
        % get the generator signal for each kernel
        cgr_generator = zeros(length(kernel_params), length(t_vals));
        for dm = 1:length(kernel_params)
            temp_kernel = Wc(:,temp_cone_indices) * kernel_params(dm).kernel';
            % make cgr generator
            for fm = 1:length(t_vals)
                mv_frame = reshape(cgr_movie(:,:,:,fm), 1, []);
                cgr_generator(dm,fm) = mv_frame * temp_kernel;
            end
            cgr_generator(dm,:) = filter(upsampled_inpulse_filter, 1, cgr_generator(dm,:));
        end
        
        
        predicted_firing_rate = multi_kernel_prediction(regressed_weights, fit_offsets, cgr_generator);
        
        figure(5)
        plot([0:pi/15:2*pi], predicted_firing_rate(34:64), 'k')
        axis([0 2*pi 0 175])
        
        L = 90;
        fz = 120; % sampling frequencey
        
        NFFT = 2^nextpow2(L);
        Y = fft(predicted_firing_rate(31:120), NFFT) / L;
        f = fz/2 * linspace(0,1,NFFT/2+1);
        power_spec = 2*abs(Y(1:NFFT/2+1));
        
        figure(11)
        plot(f(1:20), power_spec(1:20))
        
        
        f_one(sp) = sum(power_spec(4:5));
        f_two(sp) = sum(power_spec(8:9));
    end
    
    max_f_one(sf) = max(f_one);
    mean_f_one(sf) = mean(f_one);
    mean_f_two(sf) = mean(f_two);
    sd_f_one(sf) = std(f_one);
    sd_f_two(sf) = std(f_two);
    
end

freq_samples = 1./ (spatial_period *5 ./ 1000);

figure(15)
loglog(freq_samples, mean_f_two, 'r');
hold on
loglog(freq_samples, max_f_one, 'k');
hold off
title('CRG, black = F1, red = F2')
xlabel('spatial frequency, in mm')
ylabel('spike rate')
print(15,'~/Desktop/cgr.pdf', '-dpdf')

    
    
%% Drifting gratings


x_vals = 1:1:320; % stimulus width
y_vals = 1:1:320; % stimulus height
t_vals = 1:1:120; % stimulus duration

spatial_period_dg = [128 96 64 48 32 24 20 16 14 12 8 4];
%spatial_period = [96];

temporal_period = 30;

dg_f_one = zeros(length(spatial_period_dg),1);
dg_f_two = zeros(length(spatial_period_dg),1);

x = [0:6:((length(impulse_filter)-1)*6)];
xi = [0:1:((length(impulse_filter)-1)*6)];
upsampled_inpulse_filter = interp1(x, impulse_filter,xi);

for sf = 1:length(spatial_period_dg);
    
    cgr_generator = zeros(length(kernel_params), length(t_vals));
    
    
    for fm = 1:length(t_vals);
        
        stim_frame = sin((2*pi*x_vals./spatial_period_dg(sf)) + ((spatial_period_dg(sf) / temporal_period)*2*pi/spatial_period_dg(sf) * (fm-1)));
        stim_frame = repmat(stim_frame, [length(y_vals), 1, 3]);
        stim_frame = (stim_frame+1) .* 0.5;
        
        figure(3)
        image(stim_frame)
        
        mv_frame = reshape(stim_frame, 1, []);
        
        % get the generator signal for each kernel
        for dm = 1:length(kernel_params)
            temp_kernel = Wc(:,temp_cone_indices) * kernel_params(dm).kernel';
            dg_generator(dm,fm) = mv_frame * temp_kernel;
        end
    end
    
    for dm = 1:length(kernel_params)
        filtered_inputs(dm,:) = filter(upsampled_inpulse_filter, 1, dg_generator(dm,:));
    end
    
    
    predicted_firing_rate = multi_kernel_prediction(regressed_weights, fit_offsets, filtered_inputs);
    
    figure(10)
    plot(predicted_firing_rate(31:120))
    
    %[fit_amplitudes, fit_error] = fit_two_harmonics(predicted_firing_rate,...
    %                                120, temporal_period, [200 0 0]);
    
    
    L = 90;
    fz = 120; % sampling frequencey
    
    NFFT = 2^nextpow2(L);
    Y = fft(predicted_firing_rate(31:120), NFFT) / L;
    f = fz/2 * linspace(0,1,NFFT/2+1);
    
    power_spec = 2*abs(Y(1:NFFT/2+1));
    
    figure(11)
    plot(f(1:20), power_spec(1:20))
    
    
    dg_f_one(sf) = sum(power_spec(4:5));
    dg_f_two(sf) = sum(power_spec(9:10));
    
end

figure(21)
loglog(freq_samples, dg_f_one, 'k')
hold on
loglog(freq_samples, mean_f_two, 'r')
title('black = DG,F1; red = CRG, F2')
xlabel('spatial frequency, in pixels')
ylabel('spike rate')
print(21,'~/Desktop/dg_cgr.pdf', '-dpdf')

figure(22)    
loglog(freq_samples, dg_f_one, 'k') 
hold on
loglog(freq_samples, dg_f_two, 'r') 
title('black = DG,F1; red = CRG, F2')
xlabel('spatial frequency, in pixels')
ylabel('spike rate')
%print(21,'~/Desktop/dg_cgr.pdf', '-dpdf')
 
%%
x_vals = 1:1:320; % stimulus width
y_vals = 1:1:320; % stimulus height
t_vals = 1:1:120; % stimulus duration

spatial_period_dg = [128 96 64 48 32 24 20 16 14 12 8 4];
%spatial_period = [96];

temporal_period = 30;

dg_f_one_orth = zeros(length(spatial_period_dg),1);
dg_f_two_orth = zeros(length(spatial_period_dg),1);

x = [0:6:((length(impulse_filter)-1)*6)];
xi = [0:1:((length(impulse_filter)-1)*6)];
upsampled_inpulse_filter = interp1(x, impulse_filter,xi);

for sf = 1:length(spatial_period_dg);

    cgr_generator = zeros(length(kernel_params), length(t_vals));

    
    for fm = 1:length(t_vals);
        
        stim_frame = sin((2*pi*x_vals./spatial_period_dg(sf)) + ((spatial_period_dg(sf) / temporal_period)*2*pi/spatial_period_dg(sf) * (fm-1)));
        stim_frame = repmat(stim_frame, [length(y_vals), 1, 3]);
        stim_frame = (stim_frame+1) .* 0.5;
        stim_frame = permute(stim_frame, [2,1,3]);

        figure(3)
        image(stim_frame)

        mv_frame = reshape(stim_frame, 1, []);

        % get the generator signal for each kernel
        for dm = 1:length(kernel_params)
            temp_kernel = Wc(:,temp_cone_indices) * kernel_params(dm).kernel';
            dg_generator(dm,fm) = mv_frame * temp_kernel;
        end
    end
 
    for dm = 1:length(kernel_params)
        filtered_inputs(dm,:) = filter(upsampled_inpulse_filter, 1, dg_generator(dm,:));
    end


    predicted_firing_rate = multi_kernel_prediction(regressed_weights, fit_offsets, filtered_inputs);

    figure(10)
    plot(predicted_firing_rate(31:120))
    
    %[fit_amplitudes, fit_error] = fit_two_harmonics(predicted_firing_rate,...
    %                                120, temporal_period, [200 0 0]);


    L = 90;
    fz = 120; % sampling frequencey

    NFFT = 2^nextpow2(L);
    Y = fft(predicted_firing_rate(31:120), NFFT) / L;
    f = fz/2 * linspace(0,1,NFFT/2+1);

    power_spec = 2*abs(Y(1:NFFT/2+1));

    figure(11)
    plot(f(1:20), power_spec(1:20))

                                
    dg_f_one_orth(sf) = sum(power_spec(4:5));
    dg_f_two_orth(sf) = sum(power_spec(9:10));

end

figure(22)
loglog(freq_samples, (dg_f_one+dg_f_one_orth)./2, 'b')
hold on
loglog(freq_samples, mean_f_two, 'r')
title('black = DG,F1; red = CRG, F2')
xlabel('spatial frequency, in pixels')
ylabel('spike rate')
print(22,'~/Desktop/dg_cgr.pdf', '-dpdf')




%% ST-ICA


new_STE = (diag(eig_vals(1:num_sig_dims))).^0.5 * PCs(:,1:num_sig_dims)' * z_STE;

    
[ICs, A, W] = fastica(new_STE', 'g', 'gauss');


num_dims = 9;
for dm = 1:num_dims
    image_PC = Wc(:,temp_cone_indices) * PCs(:,dm);
    image_PC = reshape(image_PC, [height, width, 3]);
    image_PC = image_PC(bg_window_y:ed_window_y,bg_window_x:ed_window_x,:);
    figure(dm+2)    
    image(norm_image(image_PC))
    axis off; axis square
    title(['D',num2str(dm)])
    drawnow
end



%% clustering
num_VQs = 8;

% make response subspace from sig PCs and STA
response_subspace = [PCs(:,1:num_VQs-1), (STA'./norm(STA))];
% project stimulus into this subspace
stim_weights = STE * response_subspace;
% cluster the stimulus in this subspace

% white data along each dim
for dm = 1:size(stim_weights,2)
    stim_weights(:,dm) = stim_weights(:,dm) ./ std(stim_weights(:,dm));
end

[IDX, VQs] = kmeans(stim_weights, num_VQs, 'start', 'uniform', 'MaxIter', 500, 'OnlinePhase', 'on');

cone_space_kernels = response_subspace * VQs';

matrix_scale_factor = 1;

for comp = 1:num_VQs;
    figure(comp)
    cone_egvec = Wc(:,temp_cone_indices) * cone_space_kernels(:,comp);
    cone_egvec = reshape(cone_egvec, [320, 320, 3]);
    imagesc(matrix_scaled_up(norm_image(cone_egvec), matrix_scale_factor))
    axis off; axis image
    axis([bg_window_x ed_window_x bg_window_y ed_window_y] * matrix_scale_factor)
    title_text = ['k ', num2str(comp)];
    title(title_text)
    save_path = ['~/Desktop/cluster',num2str(comp),'.pdf'];
    %print(comp, save_path, '-dpdf', '-painters', '-r600')
    
end





