

load('/Volumes/Analysis/nishal/Demo/null_data.mat');

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
colormap gray
colorbar
%%
global_vars
stas=cell(1,1);
stas{1}=zeros(64,32,1,30);
stas{1}(32-5:32+5,16-5:16+5,1,1)=1;% rf_spatial


mov=ones(64,32,1200)*127.5;
movie_time=1200;
n_cell=1;

filt_dim1=64;
filt_dim2=32;
filt_len=30;
addpath('../lsqrSOL/');
addpath('../craigSOL/');
b= -Ax(stas,mov,movie_time,n_cell);
alpha=0.1;
beta=0.5;
momentum=0.7;

%  LSQR - Solve , change it to CRAIG when that starts working!
damp=0;
atol=10^-6;
btol=10^-6;
conlim=1.0e+300; % Doubt!
itnlim=1000;
show=1;

tic;
[ x, istop, itn, r1norm, r2norm, Anorm, Acond, Arnorm, xnorm, var ]...
   = lsqrSOL( movie_time*n_cell, movie_time*filt_dim1*filt_dim2, @dual_AAtx2, b(:), damp, atol, btol, conlim, itnlim, show );
toc;

mov_modify_new=reshape(x,[filt_dim1,filt_dim2,movie_time])+mov;
res=norm(Ax(stas,mov_modify_new,movie_time,n_cell))
mov_modify_new = mov_modify_new*127.5/max(abs(mov_modify_new(:)));
%
figure;
for itime=1200-20:1200

subplot(2,2,1);
imagesc(mov(:,:,itime));
colormap gray
colorbar
caxis([-127.5,127.5]);

subplot(2,2,2);
imagesc(mov_modify_new(:,:,itime));
colormap gray
colorbar
caxis([-127.5,127.5]);

subplot(2,2,3);
imagesc(mov(:,:,itime)-mov_modify_new(:,:,itime));
colormap gray
colorbar
caxis([-127.5,127.5]);

pause(0.1);
end

%%
global_vars
stas=cell(1,1);
stas{1}=zeros(64,32,1,30);
stas{1}(:,:,1,1)=zeros(64,32);%rf_spatial;

for itime=1:30
    stas{1}(32,16,1,itime)=(-1)^itime;
end

%mov=ones(64,32,1200)*127.5;
mov=ones(64,32,1200);
for itime=1:1200
mov(:,:,itime)=ones(64,32)*(-1)^itime*127;
end

movie_time=1200;
n_cell=1;

filt_dim1=64;
filt_dim2=32;
filt_len=30;
addpath('../lsqrSOL/');
addpath('../craigSOL/');
b= -Ax(stas,mov,movie_time,n_cell);
alpha=0.1;
beta=0.5;
momentum=0.7;

%  LSQR - Solve , change it to CRAIG when that starts working!
damp=0;
atol=10^-9;
btol=10^-9;
conlim=1.0e+300; % Doubt!
itnlim=1000;
show=1;

tic;
[ x, istop, itn, r1norm, r2norm, Anorm, Acond, Arnorm, xnorm, var ]...
   = lsqrSOL( movie_time*n_cell, movie_time*filt_dim1*filt_dim2, @dual_AAtx2, b(:), damp, atol, btol, conlim, itnlim, show );
toc;

mov_modify_new=reshape(x,[filt_dim1,filt_dim2,movie_time])+mov;
res=norm(Ax(stas,mov_modify_new,movie_time,n_cell))
mov_modify_new = mov_modify_new*127.5/max(abs(mov_modify_new(:)));
%
figure;
for itime=1:10

subplot(2,2,1);
imagesc(mov(:,:,itime));
colormap gray
colorbar
caxis([-127.5,127.5]);

subplot(2,2,2);
imagesc(mov_modify_new(:,:,itime));
colormap gray
colorbar
caxis([-127.5,127.5]);

subplot(2,2,3);
imagesc(mov(:,:,itime)-mov_modify_new(:,:,itime));
colormap gray
colorbar
caxis([-127.5,127.5]);

pause(0.1);
end

figure;
subplot(4,1,1);
ax=reshape(mov_modify_new(32,16,:),[1200,1]);
plot(ax)
subplot(4,1,2);
ax=reshape(mov_modify_new(33,16,:),[1200,1]);
plot(ax,'r')
subplot(4,1,3)
ax= reshape(mov_modify_new(32,16,:)-mov(32,16,:),[1200,1]);
plot(ax);
subplot(4,1,4)
ax= reshape(mov_modify_new(50,32,:)-mov(50,32,:),[1200,1]);
plot(ax);
