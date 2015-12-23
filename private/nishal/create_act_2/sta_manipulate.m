


for itime=4
a=full(sig_stixels);x=reshape(stas{1}(:,:,1,itime),[size(stas{1},1),size(stas{2},2)]);
subplot(2,1,1);
imagesc(x.*double(a));
caxis([min(stas{1}(:)), max(stas{1}(:))]);
colormap gray
colorbar
subplot(2,1,2);
imagesc(x);
caxis([min(stas{1}(:)), max(stas{1}(:))]);
colormap gray
colorbar
pause(1);
end

%% 

fit_params = fit_sta(stas{2});

%%
% assign the parsed results to simpler (original) names
center_point = [fit_params.center_point_x, fit_params.center_point_y];
sd_scale = [fit_params.center_sd_x, fit_params.center_sd_y];
amp_scale =1%fit_params.surround_amp_scale;
rotation_angle = fit_params.center_rotation_angle;
x_dim =fit_params.x_dim;
y_dim = fit_params.y_dim;

% initialize the output matrix
output_matrix = zeros(y_dim, x_dim);

% make an array of points for matrix (STA) values
width_points = 1:1:x_dim;
height_points = 1:1:y_dim;

% calculate the distances of these points from the center of Gauss
width_dists = center_point(1) - width_points;
height_dists = center_point(2) - height_points;

% calculate rotation matrix: couterclockwise rotation with respect to angle
rotation_matrix = [cos(rotation_angle), -1*sin(rotation_angle); sin(rotation_angle), cos(rotation_angle)];

% define covariance matrix given the sd_scale and rotation matrix
covariance_matrix = rotation_matrix * [1/sd_scale(1)^2 0; 0 1/sd_scale(2)^2] * rotation_matrix';

% calculate the value of the Gaussian at each point in output_matrix
for wd = 1:x_dim
    for ht = 1:y_dim
        pt = [height_dists(ht); width_dists(wd)];
        output_matrix(ht,wd) = amp_scale .* exp(-0.5 .* (pt' * covariance_matrix * pt));
    end
end
rf_center=output_matrix;



amp_scale =fit_params.surround_amp_scale;
surround_scale = fit_params.surround_sd_scale;
% initialize the output matrix
output_matrix = zeros(y_dim, x_dim);

% make an array of points for matrix (STA) values
width_points = 1:1:x_dim;
height_points = 1:1:y_dim;

% calculate the distances of these points from the center of Gauss
width_dists = center_point(1) - width_points;
height_dists = center_point(2) - height_points;

% calculate rotation matrix: couterclockwise rotation with respect to angle
rotation_matrix = [cos(rotation_angle), -1*sin(rotation_angle); sin(rotation_angle), cos(rotation_angle)];

% define covariance matrix given the sd_scale and rotation matrix
covariance_matrix = rotation_matrix * [1/sd_scale(1)^2 0; 0 1/sd_scale(2)^2] * rotation_matrix';

% calculate the value of the Gaussian at each point in output_matrix
for wd = 1:x_dim
    for ht = 1:y_dim
        pt = surround_scale*[height_dists(ht); width_dists(wd)]; %% Nishal DOUBT TODO
        output_matrix(ht,wd) = amp_scale .* exp(-0.5 .* (pt' * covariance_matrix * pt));
    end
end

rf_surround=output_matrix;

rf_spatial=rf_center-rf_surround;

scale_one=fit_params.scale_one;
scale_two=fit_params.scale_two;
tau_one=fit_params.tau_one;
tau_two=fit_params.tau_two;
n_filters=fit_params.n_filters;
t=[0:29];
tf = scale_one*((t/tau_one).^n_filters).*exp(-n_filters*(t/tau_one -1)) - scale_two*((t/tau_two).^n_filters).*exp(-n_filters*(t/tau_two -1));
figure;plot(tf)

figure;
imagesc(rf_spatial);


%%

datafile = '2013-10-10-0/data000';
type_name= cell(1,1);
type_name{1}='On Parasol';

datarun=load_data(datafile)
datarun=load_sta(datarun)
datarun=load_params(datarun)

%get_cell_ids(datarun,type_name) % Vision IDs - used for loading fitted
%STAs!

matlab_cell_ids=get_cell_indices(datarun,type_name);
stas=datarun.stas.stas(matlab_cell_ids);
n_cell=length(stas);

datarun = compute_sta_fits(datarun, type_name{1})

addpath('../fwdfittingfunctions/');

sta_fit_cs=cell(21,1);

for icell=1:66
    icell
sta_fit_cs{icell}  = fit_sta_sequence(stas{icell},  'fit_temporal',false,'fit_center',true,'fit_surround',true);
end















