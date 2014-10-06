% this script is meant for testing during the development of fit_sta

%datarun = load_data('2008-08-27-6', 'rf-9-js-1')
datarun = load_data('2008-08-27-5', 'rf-8-gf');
%datarun = load_data('/snle/lab/Experiments/Array/Analysis/2008-08-27-6/data009/data009');
%datarun = load_data('/snle/lab/Experiments/Array/Analysis/2008-08-27-6/data017/data017');
datarun = load_data('2009-04-13-5', 'rf-8-gf-center-obvius-stas');


datarun = load_data('2009-04-13-5', 'rf-8-gf-center-obvius-stas');

datarun = load_data('2009-04-13-5', 'rf-8-gf-center-obvius-stas');

datarun = load_data('2008-08-27-6', 'rf-0-gf');
datarun = load_data('2008-08-27-6', 'rf-0-gf-vision-stas');


datarun = load_params(datarun);
datarun = load_sta(datarun);
datarun = load_neurons(datarun);

datarun = load_obvius_sta_fits(datarun);
datarun = get_sta_summaries(datarun, {1,2,3,4});
datarun = get_sta_fits_from_obvius(datarun, 'all');

%% pick a particular cell out 

cell_id = 6;
plot_rf(datarun, cell_id)
temp_index = get_cell_indices(datarun, cell_id);


sta = datarun.stas.stas{temp_index};
tmp_frames = size(sta,4);
if datarun.stimulus.interval == 1
    sta = sta(:,:,:,1:(tmp_frames-1));
end
rf = datarun.stas.rfs{temp_index};


%% make a fake STA

fake_params = [14.5,...
                18.25,...
                3.25,...
                0.65,...
                -0.5,...
                0.3,...
                1,...
                0.33,...
                64,...
                32,...
                3,...
                0,...
                0,...
                0,...
                -0.2,...
                0.05,...
                3.5,...
                11.25,...
                6,...
                30];     
            
fake_sta = sta_fit_function(fake_params);

% get peak value of fake_sta
peak_sta = max(abs(reshape(fake_sta, 1, [])));
% make noise that is a few percent of the peak sta value (mimics a high SNR
% sta
sta_noise = normrnd(0, peak_sta * 0.02, [32, 64, 3, 30]);
fake_sta = reshape(fake_sta, [],1) + reshape(sta_noise,[],1);
fake_sta = reshape(fake_sta, [32, 64, 3, 30]);

fake_sig_stixels = significant_stixels(fake_sta, 'time', 'max');
fake_rf = rf_from_sta(fake_sta, 'sig_stixels', fake_sig_stixels);
fake_time_course = time_course_from_sta(fake_sta, fake_sig_stixels);

figure(20)
imagesc(norm_image(fake_rf))

figure(21)
plot(fake_time_course(:,1), 'r')
hold on
plot(fake_time_course(:,2), 'g')
plot(fake_time_course(:,3), 'b')
hold off

%%

[fit_params, initial_params, sta_fit] = fit_sta(fake_sta, 'fit_surround', false,...
                'fit_surround_amp_scale', false,...
                'fit_surround_sd_scale', false,...
                'initial_center_point_x', 14.5,...
                'initial_center_point_y', 18.25,...
                'initial_tau_one', 3.5,...
                'initial_tau_two', 11.25,...
                'initial_scale_one', -0.2,...
                'initial_scale_two', 0.05,...
                'initial_color_weight_a', 0.3,...
                'initial_color_weight_b', 1,...
                'initial_color_weight_c', 0.33,...
                'fit_n_filters', false,...
                'fit_color_weight_a', false,...
                'fit_color_weight_b', false,...
                'fit_color_weight_c', false,...
                'fit_scale_one', false,...
                'fit_scale_two', false,...
                'fit_tau_two', true,...
                'fit_tau_one', true,...
                'fit_center_point_x', true,...
                'fit_center_point_y', true,...
                'initial_n_filters', 6);

            
[fit_ = fit_sta(fake_sta, 'fit_surround', false)

fit_params

param_mat = zeros(24,20);
tic
    for cnt = 1:8
        
    param_mat(cnt,:) = fit_sta(fake_sta, 'fit_surround', false,...
                    'fit_surround_amp_scale', false,...
                    'fit_surround_sd_scale', false,...
                    'fit_n_filters', false,...
                    'initial_n_filters', 6);
    %param_mat(cnt,:) = fit_params;  
    
    end
toc



% FIT INFORMATION 
% 1 center_point_x
% 2 center_point_y
% 3 sd_x
% 4 sd_y
% 5 rotation_angle
% 6 color_weight_a
% 7 color_weight_b
% 8 color weight_c
% 9 x_dim
% 10 y_dim
% 11 num_colors
% 12 fit_surround
% 13 surround_sd_scale
% 14 surround_amp_scale (relative to center_amp)
% 15 scale_one
% 16 scale_two
% 17 tau_one
% 18 tau_two
% 19 n-filters
% 20 frame_number
    
%%

% view the rf of interest

mark_selection.thresh = 4.0;

fit_params = fit_sta(sta, 'fit_surround', false,...
                'fit_surround_amp_scale', false,...
                'fit_surround_sd_scale', false,...
                'fit_n_filters', false,...
                'verbose', false,...
                'mark_params', mark_selection)
            

            
fit_params = fit_sta(sta, 'fit_surround', true,...
                'fit_surround_amp_scale', true,...
                'fit_surround_sd_scale', false,...
                'initial_tau_one', [],... 
                'initial_tau_two', [],... 
                'initial_scale_one', [],...
                'initial_scale_two', [],...
                'initial_color_weight_a', datarun.obvius.sta_fits{temp_index}.color(1),...
                'initial_color_weight_b', datarun.obvius.sta_fits{temp_index}.color(2),...
                'initial_color_weight_c', datarun.obvius.sta_fits{temp_index}.color(3),...
                'initial_center_point_x', datarun.stas.fits{temp_index}.mean(1),...
                'initial_center_point_y', datarun.stas.fits{temp_index}.mean(2),...
                'initial_surround_amp_scale', -1*datarun.stas.fits{temp_index}.surround_scale,...
                'initial_surround_sd_scale', 2.0,...
                'initial_center_sd_x', datarun.stas.fits{temp_index}.sd(2),...
                'initial_center_sd_y', datarun.stas.fits{temp_index}.sd(1),...
                'initial_center_rotation_angle', datarun.stas.fits{temp_index}.angle,...
                'fit_center_sd_x', true,...
                'fit_center_sd_y', true,...
                'fit_center_rotation_angle', true,...
                'fit_scale_one', true,...
                'fit_scale_two', true,...
                'fit_tau_two', true,...
                'fit_tau_one', true,...
                'fit_color_weight_a', true,...
                'fit_color_weight_b', true,...
                'fit_color_weight_c', true,...
                'fit_center_point_x', true,...
                'fit_center_point_y', true,...
                'initial_n_filters', 6,...
                'fit_n_filters', false,...
                'verbose', false);
              

fit_params
initial_params

plot_rf_summaries(datarun, cell_id, 'foa', 2, 'clear', false, 'plot_fits', true, 'fit_color', 'r');



%sta_fit = make_fit_sta(fit_params);
temp_sig_stixels = significant_stixels(sta_fit, 'time', 'max', 'thresh', 3);
real_sig_stixels = significant_stixels(sta, 'time', 'max', 'thresh', 3);
fit_time_course = time_course_from_sta(sta_fit, temp_sig_stixels);
figure(12)
plot(fit_time_course ./ abs(ext(fit_time_course(:,2))), '--', 'LineWidth', 2)
hold on
temp_tcs = datarun.stas.time_courses{1};
plot(temp_tcs ./ abs(ext(temp_tcs(:,2))))
hold off


rf_fit = rf_from_sta(sta_fit, 'sig_stixels', temp_sig_stixels);
rf_fit = rf_fit ./ abs(ext(reshape(rf_fit, 1, [])));
figure(2)
imagesc(norm_image(rf_fit)) 
axis image

figure(1)
rf = rf_from_sta(sta, 'sig_stixels', real_sig_stixels);
rf = rf ./ abs(ext(reshape(rf, 1, [])));
imagesc(norm_image(rf))
axis image

residual_rf = rf - rf_fit;
figure(3)
imagesc(norm_image(residual_rf))
axis image

figure(5)
imagesc(norm_image(sta_fit(:,:,:,25)))


%%

sig_stix(1:32,1:64) = false;
sig_stix(14,25) = true;
sta_tc = time_course_from_sta(sta,sig_stix);
figure(1)
hold on
plot(sta_tc(:,1), 'r')
plot(sta_tc(:,2), 'g')
plot(sta_tc(:,3), 'b')





%% figure out how time course fit works

t_points = 0:1:50;
p_one = -0.20;
p_two = 0.054;
tau_one = 3.90;
tau_two = 11.17;
n = 6;

filt_one = p_one .* ((t_points ./ tau_one).^n) .* exp((-1 .* n) .* ((t_points ./tau_one) - 1));
filt_two = p_two .* ((t_points ./ tau_two).^n) .* exp((-1 .* n) .* ((t_points ./tau_two) - 1));
tc = filt_one + filt_two;

figure(4)
plot(t_points, tc)


%%
figure
tmp_index = get_cell_indices(datarun, 1);
tmp_tc = datarun.stas.time_courses{tmp_index};
plot(tmp_tc(:,1), 'r')
hold on
plot(tmp_tc(:,2), 'g')
plot(tmp_tc(:,3), 'b')
hold off

%%
figure(10);clf
drawEllipse([[0 0] [0.75 2*0.9] (-0.813-pi/2)])
axis([-10 10 -10 10])


    
%%    

cell_type = {1};

clear fit_params

mark_selection.thresh = 5.0;

fit_params.fit_surround = false;
fit_params.fit_surround_sd_scale = false;
fit_params.fit_surround_amp_scale = false;
fit_params.mark_params = mark_selection;


% If stimulus interval = 1, shave off a frame to correct for bug in vision
if datarun.stimulus.interval == 1
    corrected_datarun = datarun;
    num_frames = size(datarun.stas.stas{1}, 4);
    for rgc = 1:length(datarun.cell_ids);
        temp_sta = datarun.stas.stas{rgc};
        temp_sta = temp_sta(:,:,:,1:(num_frames -1));
        corrected_datarun.stas.stas{rgc} = temp_sta;
    end
else
    corrected_datarun = datarun;
end
            

corrected_datarun = compute_sta_fits_multicore(corrected_datarun, cell_type, 'fit_params', fit_params);


cell_indices = get_cell_indices(corrected_datarun,cell_type);
new_fit_datarun = corrected_datarun;
fit_struct = new_fit_datarun.matlab.sta_fits;
for cc = 1:length(cell_indices)
    new_fit_datarun.stas.fits{cell_indices(cc)}.mean = [fit_struct{cell_indices(cc)}.center_point_x,...
                                                fit_struct{cell_indices(cc)}.center_point_y];
                                            
    new_fit_datarun.stas.fits{cell_indices(cc)}.sd = [fit_struct{cell_indices(cc)}.center_sd_y,...
                                              fit_struct{cell_indices(cc)}.center_sd_x];
                                          
    new_fit_datarun.stas.fits{cell_indices(cc)}.angle = -fit_struct{cell_indices(cc)}.center_rotation_angle;
    new_fit_datarun.stas.fits{cell_indices(cc)}.surround_sd_scale = fit_struct{cell_indices(cc)}.surround_sd_scale;
    new_fit_datarun.stas.fits{cell_indices(cc)}.surround_scale = fit_struct{cell_indices(cc)}.surround_amp_scale;;
end

%corrected_datarun = compute_sta_fits_multicore(corrected_datarun, 6, 'fit_params', fit_params);


    
%% compare spatial fits between obvius and matlab
cell_type = {1};

plot_rf_summaries(datarun, cell_type, 'foa', 5, 'clear', true,...
                     'plot_fits', true, 'fit_color', 'r', 'label', true);

plot_rf_summaries(new_fit_datarun, cell_type, 'foa', 5, 'clear', false,...
                        'plot_fits', true, 'fit_color', 'k', 'label', true);

%vision_datarun = datarun;
%vision_datarun = get_sta_fits_from_vision(vision_datarun);
plot_rf_summaries(vision_datarun, cell_type, 'foa', 5, 'clear', false,...
                        'plot_fits', true, 'fit_color', 'k', 'label', true);
    
    
%% compare temporal fits between obvius and matlab

cell_type = {1};

temp_indices = get_cell_indices(new_fit_datarun, cell_type);
num_rgcs = length(temp_indices);


tau_one_comp = zeros(num_rgcs,1);
tau_two_comp = zeros(num_rgcs,1);
scale_one_comp = zeros(num_rgcs,1);
scale_two_comp = zeros(num_rgcs,1);
for cc = 1:num_rgcs
    tau_one_obvius = new_fit_datarun.obvius.sta_fits{temp_indices(cc)}.tau(1);
    tau_two_obvius = new_fit_datarun.obvius.sta_fits{temp_indices(cc)}.tau(2);
    
    scale_one_obvius = new_fit_datarun.obvius.sta_fits{temp_indices(cc)}.scale(1);
    scale_two_obvius = new_fit_datarun.obvius.sta_fits{temp_indices(cc)}.scale(2);

    tau_one_matlab = new_fit_datarun.matlab.sta_fits{temp_indices(cc)}.tau_one;
    tau_two_matlab = new_fit_datarun.matlab.sta_fits{temp_indices(cc)}.tau_two;
    
    scale_one_matlab = new_fit_datarun.matlab.sta_fits{temp_indices(cc)}.scale_one;
    scale_two_matlab = new_fit_datarun.matlab.sta_fits{temp_indices(cc)}.scale_two;
    
    tau_one_comp(cc) = tau_one_matlab ./ tau_one_obvius;
    tau_two_comp(cc) = tau_two_matlab ./ tau_two_obvius;
    scale_one_comp(cc) = scale_one_matlab ./ scale_one_obvius;
    scale_two_comp(cc) = scale_two_matlab ./ scale_two_obvius;
end

mean(tau_one_comp)
mean(tau_two_comp)
mean(scale_one_comp)
mean(scale_two_comp)


figure(33)
subplot(2,2,1)
[tau_one_hist, tau_one_bins] = hist(tau_one_comp);
bar(tau_one_bins, tau_one_hist)
xlabel('tau one')

subplot(2,2,2)
[tau_two_hist, tau_two_bins] = hist(tau_two_comp);
bar(tau_two_bins, tau_two_hist)
xlabel('tau two')

subplot(2,2,3)
[scale_one_hist, scale_one_bins] = hist(scale_one_comp);
bar(scale_one_bins, scale_one_hist)
xlabel('scale one')

subplot(2,2,4)
[scale_two_hist, scale_two_bins] = hist(scale_two_comp);
bar(scale_two_bins, scale_two_hist)    
xlabel('scale two')
    
    
fit_time_course = sta_fit_time_course(new_fit_datarun.matlab.sta_fits{3}, 'matlab');
figure(35); clf; hold on
plot(fit_time_course(:,1), 'r')
plot(fit_time_course(:,2), 'g')
plot(fit_time_course(:,3), 'b')


fit_time_course = sta_fit_time_course(new_fit_datarun.obvius.sta_fits{3}, 'obvius');
plot(fit_time_course(:,1), '--r')
plot(fit_time_course(:,2), '--g')
plot(fit_time_course(:,3), '--b')

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    








