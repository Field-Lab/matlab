% Try different fits and see performance .. 
%%
startup_bertha
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

% stas manipulate ? 
stas_new=cell(length(stas),1);
for icell=1:length(stas)
    st_temp=zeros(size(stas{1},2),size(stas{1},1),1,size(stas{1},4)); % DOUBT .. Could be a reason for things to fail!!!!!
    for itime=1:30
        st_temp(:,:,:,itime)=mean(stas{icell}(:,:,:,end-itime+1),3)'; % DOUBT .. Could be a reason for things to fail!!!!!
    end
    stas_new{icell}=st_temp;
end
stas=stas_new;
clear stas_new
filt_len=size(stas{1},4);

%
filt_dim1=size(stas{1},1);
filt_dim2=size(stas{1},2);

stas_fit=cell(n_cell,1);

addpath('../fwdfittingfunctions/');
% Add matlab code path addpath(genpath('../../../code/'));
datarun = compute_sta_fits(datarun, type_name{1});

%% Extract approximated STA from datarun

% assign the parsed results to simpler (original) names
for icell=1:n_cell
fit_params=datarun.matlab.sta_fits{matlab_cell_ids(icell)};

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
rf_spatial=rf_spatial';
%rf_spatial=rf_spatial(end:-1:1,end:-1:1);

scale_one=fit_params.scale_one;
scale_two=fit_params.scale_two;
tau_one=fit_params.tau_one;
tau_two=fit_params.tau_two;
n_filters=fit_params.n_filters;
t=[0:29];
tf = scale_one*((t/tau_one).^n_filters).*exp(-n_filters*(t/tau_one -1)) - scale_two*((t/tau_two).^n_filters).*exp(-n_filters*(t/tau_two -1));


% do flipping of RF ? 
sta_fitted = zeros(filt_dim1,filt_dim2,1,filt_len);
for itime=1:filt_len
sta_fitted(:,:,1,itime)= rf_spatial * tf(itime);
end

stas_fit{icell} = sta_fitted;


subplot(2,2,1);
imagesc(squeeze(stas{icell}(:,:,1,3)));
colormap gray
colorbar
axis image

subplot(2,2,2);
imagesc(rf_spatial*max(tf));
colormap gray
colorbar
axis image

subplot(2,2,3);
plot(tf)

pause(1)

end

% Play movie for an STA
figure;
icell=3;
for itime=1:filt_len
    itime
    subplot(2,1,1);
    imagesc(stas{icell}(:,:,1,itime));
     caxis([min(stas{icell}(:)) max(stas{icell}(:))]);
    colormap gray
    colorbar
    axis image
   
    subplot(2,1,2);
    imagesc(stas_fit{icell}(:,:,1,itime));
    caxis([min(stas{icell}(:)) max(stas{icell}(:))]);
    colormap gray
    colorbar
    axis image

    
    pause(0.5)
end

%% Just Thereshold STAs 
stas_thr=cell(n_cell,1);

for icell=1:n_cell
    % clip outside 5 sigma 
    fit_params=datarun.matlab.sta_fits{matlab_cell_ids(icell)};

center_point = [fit_params.center_point_x, fit_params.center_point_y];
sd_scale = [fit_params.center_sd_x, fit_params.center_sd_y];
amp_scale =1%fit_params.surround_amp_scale;
rotation_angle = fit_params.center_rotation_angle;
x_dim =fit_params.x_dim;
y_dim = fit_params.y_dim;

% initialize the output matrix
dist_matrix = zeros(y_dim, x_dim);

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
        dist_matrix(ht,wd) =  sqrt(pt' * covariance_matrix * pt);
    end
end

% 
% 
% amp_scale =fit_params.surround_amp_scale;
% surround_scale = fit_params.surround_sd_scale;
% % initialize the output matrix
% output_matrix = zeros(y_dim, x_dim);
% 
% % make an array of points for matrix (STA) values
% width_points = 1:1:x_dim;
% height_points = 1:1:y_dim;
% 
% % calculate the distances of these points from the center of Gauss
% width_dists = center_point(1) - width_points;
% height_dists = center_point(2) - height_points;
% 
% % calculate rotation matrix: couterclockwise rotation with respect to angle
% rotation_matrix = [cos(rotation_angle), -1*sin(rotation_angle); sin(rotation_angle), cos(rotation_angle)];
% 
% % define covariance matrix given the sd_scale and rotation matrix
% covariance_matrix = rotation_matrix * [1/sd_scale(1)^2 0; 0 1/sd_scale(2)^2] * rotation_matrix';
% 
% % calculate the value of the Gaussian at each point in output_matrix
% for wd = 1:x_dim
%     for ht = 1:y_dim
%         pt = surround_scale*[height_dists(ht); width_dists(wd)]; %% Nishal DOUBT TODO
%         output_matrix(ht,wd) = amp_scale .* exp(-0.5 .* (pt' * covariance_matrix * pt));
%     end
% end
% 
% rf_surround=output_matrix;
% 
% rf_spatial=rf_center-rf_surround;
% rf_spatial=rf_spatial';
% %rf_spatial=rf_spatial(end:-1:1,end:-1:1);

dist_matrix=dist_matrix';
% do flipping of RF ? 
sta_th = zeros(filt_dim1,filt_dim2,1,filt_len);
for itime=1:filt_len
sta_th(:,:,1,itime)= stas{icell}(:,:,1,itime).*(dist_matrix<=5);% as 5 sigma
end

stas_thr{icell} = sta_th;

figure;
subplot(1,2,1);
imagesc(squeeze(stas{icell}(:,:,1,3)));
colormap gray
colorbar
axis image

subplot(1,2,2);
imagesc(squeeze(stas_thr{icell}(:,:,1,3)));
colormap gray
colorbar
axis image

pause(1)


end

% Have stas, stas_thr,stas_fit
%%
[mov_orig1,mov_modify_new1] = compute_for_stas(stas);
[mov_orig2,mov_modify_new2] = compute_for_stas(stas_fit);
[mov_orig3,mov_modify_new3] = compute_for_stas(stas_thr);

%% 

cell_resp11=Ax(stas,mov_modify_new1,size(mov_modify_new1,3),n_cell);
cell_resp21=Ax(stas_fit,mov_modify_new1,size(mov_modify_new1,3),n_cell);
cell_resp31=Ax(stas_thr,mov_modify_new1,size(mov_modify_new1,3),n_cell);

figure;
subplot(3,1,1);
plot(cell_resp11);
subplot(3,1,2);
plot(cell_resp21);
subplot(3,1,3);
plot(cell_resp31);


cell_resp12=Ax(stas,mov_modify_new2,size(mov_modify_new2,3),n_cell);
cell_resp22=Ax(stas_fit,mov_modify_new2,size(mov_modify_new2,3),n_cell);
cell_resp32=Ax(stas_thr,mov_modify_new2,size(mov_modify_new2,3),n_cell);

figure;
subplot(3,1,1);
plot(cell_resp12);
subplot(3,1,2);
plot(cell_resp22);
subplot(3,1,3);
plot(cell_resp32);

cell_resp13=Ax(stas,mov_modify_new3,size(mov_modify_new3,3),n_cell);
cell_resp23=Ax(stas_fit,mov_modify_new3,size(mov_modify_new3,3),n_cell);
cell_resp33=Ax(stas_thr,mov_modify_new3,size(mov_modify_new3,3),n_cell);

figure;
subplot(3,1,1);
plot(cell_resp13);
subplot(3,1,2);
plot(cell_resp23);
subplot(3,1,3);
plot(cell_resp33);




