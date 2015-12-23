load('/Volumes/Analysis/nishal/data/sta_fit_white_null_noise.mat');
load('/Volumes/Analysis/nishal/data/white_null_mov2.mat','stas')

% Have datarun now!
%%

A_mat=[];
figure;
icnt=0;
for icell=1:length(datarun.matlab.sta_fits)
if(~isempty(datarun.matlab.sta_fits{icell}))
        icnt=icnt+1;
fit_params= sta_fit_cs{icnt};% datarun.matlab.sta_fits{icell};
    icell

    % assign the parsed results to simpler (original) names
center_point = [fit_params.center_point_x, fit_params.center_point_y];
sd_scale = [fit_params.center_sd_x, fit_params.center_sd_y];
amp_scale =1
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


subplot(2,1,1);
imagesc(stas{icnt}(:,:,1,28));
colormap gray
colorbar

subplot(2,1,2);
imagesc(rf_center);
%title(sprintf('Vision ID %d',datarun.cell_ids(icell)));
colormap gray
colorbar
pause(0.5);
A_mat=[A_mat;rf_spatial(:)'];
%x=stas{icnt}(:,:,1,5);
%A_mat=[A_mat;x(:)'];
end

end

%%

var64=64;
filt_dim1=var64;
filt_dim2=32;
movie_white_len=120*100;
mov=zeros(filt_dim1,filt_dim2,movie_white_len);
init_mov=zeros(filt_dim1,filt_dim2,movie_white_len);
n_cell=size(A_mat,1);
filt_len=30;

tic;
A_right_inv = (A_mat'*((A_mat*A_mat')\eye(n_cell,n_cell)));
toc;

tic;
for itime=filt_len:movie_white_len

 
    
    init_frame= randn(filt_dim1,filt_dim2);
    b = - A_mat*init_frame(:) ;
    x=A_right_inv*b;
    
% LSQR could be used here .. still have to implement Atx     
% damp=0;
% atol=10^-6;
% btol=10^-6;
% conlim=1.0e+300; % Doubt!
% itnlim=1000;
% show=1;

% tic;
% [ x, istop, itn, r1norm, r2norm, Anorm, Acond, Arnorm, xnorm, var1 ]...
%    = lsqrSOL( n_cell, filt_dim1*filt_dim2, @dual_AAtx2_frame, b(:), damp, atol, btol, conlim, itnlim, show );
% toc;
% Or, do simple matrix vector multiplication! .. and send to GPU in that
% case!


mov(:,:,itime)=reshape(x,[filt_dim1,filt_dim2])+init_frame;
init_mov(:,:,itime)=init_frame;
end
toc;

%%
figure;
for itime=1:120:movie_white_len
    subplot(2,1,1);
    imagesc(init_mov(:,:,itime));
    colormap gray
    colorbar
    
    subplot(2,1,2);
    imagesc(mov(:,:,itime));
    colormap gray
    colorbar
    %caxis([-1,2]);
    %title(sprintf('Time: %d',itime));
    pause(0.1);
end
%%
tlen=12000;
figure;
subplot(2,1,1)
nonnull=(Ax(stas,init_mov(:,:,1:tlen),tlen,n_cell));
plot(nonnull);
title('White movie')
subplot(2,1,2)
nulll=(Ax(stas,mov(:,:,1:tlen),tlen,n_cell));
plot(nulll);
title('Spatial nulling')

cellR=zeros(n_cell,tlen);
for itime=1:tlen
    x=mov(:,:,itime);
cellR(:,itime)=A_mat*x(:);

end
figure;
plot(cellR')


