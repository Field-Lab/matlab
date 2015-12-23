%%
data_paths{1} = '2007-03-02-1/data016-gdf/data016/data016';  % cone
data_paths{1} = '2007-03-02-1/data001-map-gdf/data001-map';  % rod




%%
num_sets = length(data_paths);
cell_types = {2};

for dset = 1:num_sets
    datarun{dset} = load_data(data_paths{dset});
    datarun{dset} = load_sta(datarun{dset}, 'load_sta', []);
    datarun{dset} = load_params(datarun{dset});
    datarun{dset} = load_neurons(datarun{dset});

    marks_params.thresh = 4.0;  % set the threshold for significant pixels (units sd)
    datarun{dset} = get_sta_summaries(datarun{dset}, cell_types, 'marks_params', marks_params);
end


%plot_average_rf_profile(datarun{1},{2})

% compute the average RF
scale = 3;
ave_rf = get_average_rf(datarun{1}, {2}, 'scale', scale);

% get the sig stix and com
sig_stixels = significant_stixels(ave_rf);
ave_rf_com = rf_com(ave_rf, 'sig_stixels', sig_stixels);
image_center = round(ave_rf_com);

window_radius = 15; % rod
window_radius = 35; % cone

x_image_begin = image_center(2) - window_radius;
x_image_end = image_center(2) + window_radius;
y_image_begin = image_center(1) - window_radius;
y_image_end = image_center(1) + window_radius;
windowed_rf = ave_rf(x_image_begin:x_image_end,y_image_begin:y_image_end,:); 

figure(5)
imagesc(norm_image(-1*matrix_scaled_up(windowed_rf, 10)))
axis equal
print(5,'~/Desktop/spatial_rf.pdf', '-dpdf')


% figure(1); clf;
% plot(temp_profile(:,:,1), 'k')
% hold on
% plot(temp_profile(:,:,2), 'g')
% plot(temp_profile(:,:,3), 'b')
% axis([36 80 -0.1 1])
% print(1, '~/Desktop/photopic-profiles', '-dpdf')

windowed_sig_stixels = significant_stixels(windowed_rf);
windowed_rf_com = rf_com(windowed_rf, 'sig_stixels', windowed_sig_stixels);

if 0
    [ave_profile_x, ave_profile_y] = rf_profile(windowed_rf, windowed_rf_com);
else
    ave_profile_y = squeeze(windowed_rf(30,:,:));
    %ave_profile_y = squeeze(windowed_rf(round(windowed_rf_com(1)),:,:));
    ave_profile_x = 1:1:max(size(ave_profile_y));
end

figure(1)
if size(ave_profile_y,2) == 3
    plot(ave_profile_x, ave_profile_y(:,1) ./ max(ave_profile_y(:,2)),'r')
    hold on
    plot(ave_profile_x, ave_profile_y(:,2) ./ max(ave_profile_y(:,2)),'g')
    plot(ave_profile_x, ave_profile_y(:,3) ./ max(ave_profile_y(:,2)),'b')
    hold off
    axis([1 length(ave_profile_x) -0.4 1.4])
    print(1, '~/Desktop/color_profile.pdf', '-dpdf')
else
    plot(ave_profile_x, ave_profile_y ./ max(ave_profile_y), 'k')
    axis([1 length(ave_profile_x) -0.4 1.4])
    print(1, '~/Desktop/color_profile.pdf', '-dpdf')
end

if size(ave_profile_y,2) == 3
    ave_profile_y_fit = ave_profile_y(:,2) ./ max(ave_profile_y(:,2));
else
    ave_profile_y_fit = ave_profile_y ./ max(ave_profile_y);
end

fit_params.fit_surround_radius = true;
fit_params.fit_surround_scale = true;
fit_params.fit_center = true;
fit_params.center = 30;
fit_params.center_scale = 1;
fit_params.center_radius = 5;
fit_params.surround_radius = 10;
fit_params.surround_scale = 2;
[fit_params, the_fit] = fit_profile(ave_profile_x,ave_profile_y_fit, fit_params);
fit_params

figure(2)
plot(ave_profile_x, ave_profile_y_fit)
hold on
plot(the_fit(:,1), the_fit(:,2), 'r')
hold off
axis([1 length(ave_profile_x) -0.4 1.4])
print(2,'~/Desktop/profile_and_fit.pdf', '-dpdf')


center_gauss = normpdf(ave_profile_x,fit_params.center,fit_params.center_radius);
center_gauss = fit_params.center_scale * (center_gauss ./ max(center_gauss));
surround_gauss = normpdf(ave_profile_x,fit_params.center, fit_params.surround_radius);
surround_gauss = -1*fit_params.surround_scale * (surround_gauss ./ max(surround_gauss));
figure(3)
plot(ave_profile_x, center_gauss, 'r', ave_profile_x, surround_gauss, 'k')
axis([1 length(ave_profile_x) -0.4 1.4])
print(3, '~/Desktop/center_surround.pdf', '-dpdf')



% try to fit cone profile with rod center

fit_params.fit_surround_radius = true;
fit_params.fit_surround_scale = true;
fit_params.fit_center = true;
fit_params.fit_center_radius = false;
fit_params.center = 36;
fit_params.center_scale = 4;
fit_params.center_radius = 3.4330 * 2
fit_params.surround_radius = 30;
fit_params.surround_scale = 3;
[fit_params, the_fit] = fit_profile(ave_profile_x,ave_profile_y_fit, fit_params);
fit_params

figure(2)
plot(ave_profile_x, ave_profile_y_fit)
hold on
plot(the_fit(:,1), the_fit(:,2), 'r')
hold off
axis([1 length(ave_profile_x) -0.4 1.4])
print(2,'~/Desktop/profile_and_fit.pdf', '-dpdf')






