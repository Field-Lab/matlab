mov_gen2=zeros(Filtdim1,Filtdim2,movieLen);
movie_idx=2;


if(movie_idx==1)
mov_gen2(319-310,158-150,:)=0.5; % Stimulate a sub-unit
mov_gen2(320-310,160-150,:)=0.5;
mov_gen2=mov_gen2+rand(Filtdim1,Filtdim2,movieLen)*0.3;
end

if(movie_idx==2)
mov_gen2=double(rand(Filtdim1,Filtdim2,movieLen)>0.5)-0.5;
end
if(movie_idx==3)
    latency=30;
mov_short=double(rand(Filtdim1,Filtdim2,ceil(movieLen/latency)+2)>0.5)-0.5;
mov_gen2=zeros(Filtdim1,Filtdim2,movieLen);
    for iframe=1:movieLen
        mov_gen2(:,:,iframe)=mov_short(:,:,floor(iframe/latency)+1);
    end
end
if(movie_idx==4)
    latency=20;
mov_short=double(rand(Filtdim1,Filtdim2,ceil(movieLen/latency)+2));
mov_gen2=zeros(Filtdim1,Filtdim2,movieLen);
    for iframe=1:movieLen
        mov_gen2(:,:,iframe)=mov_short(:,:,floor(iframe/latency)+1);
    end
end

mov_gen2(:,:,1:30)=0;

figure;
for itime=40:50
    itime
imagesc(mov_gen2(:,:,itime));
colormap gray
colorbar
axis image
caxis([-0.5,0.5]);
pause(0.01)

end


addpath(genpath('../create_act_2/'));

mov2=zeros(Filtdim1,Filtdim2,movieLen+2*120);
mov2(:,:,121:end-120)=mov_gen2;

mov=(127.5/0.5)*mov2;clear mov2
mov_orig2 = mov;
%% 
cell_params=struct();
cell_params.type_name_inp='OFF Parasol';%'userCellList';
cell_params.cell_list=[];
% cell_params.STAlen=14;
% cell_params.sta_spatial=sprintf('%s/stas_spatial.mat',destination_mat);
cell_params.use_fits=2; % 2, 0,0,2
cell_params.sta_spatial_method=4;%1,2 ,3,4
% cell_params.sta_spatial_method = 1 for just using 4th frame, 2 is for fitting spatial STA. 
% Use cell_params.use_fits=2 (clipping) if cell_params.sta_spatial_method = 1 and 
% use cell_params.use_fits=0 (no processing of STA) if
% cell_params.sta_spatial_method = 2;
% STA spatial null Method 3 = low rank, 4 = average waveform and use it ..  

mov_params=struct();
mov_params.mov_type='bw-precomputed';
mov_params.movie_time=120*15*60/2;
mov_params.mean=0.5*255;
mov_params.deviation=0.48*255;

mov_params.mdf_file = '/Volumes/Analysis/stimuli/white-noise-xml/BW-8-2-0.48-11111.xml';
mov_params.stixel=8;

mov_params.interval = 2; % Not important to have this parameter. Default is 1. When we want repeated frames, just set this interval (This just controls the blank Frames), and select the movie_time appropriately.

% Post process. Default is stretch. If using default, need to give only mov_params.scaling_loss parameter.
mov_params.post_process_method = 'scale'; % or, 'stretch'
mov_params.scale = 0.48/0.48;
%mov_params.scaling_loss=0.05; % a number in [0,1], fraction of values that is changed by scaling.

solver=21; % Solver 4 used for spatial nulling, 7 for iterated spatial nulling


matlab_cell_ids=0;
cell_params.CellMasks = CellMasks;
if(solver==21) % Dykstra's Alternating Projections, Spatial, GPU - be on smokestack
     togo=1;
   %  mov=gpuArray(mov);mov_orig = gpuArray(mov_orig);
    while togo==1
    maxClip = (0.48/0.5)*127.5;
[~,mov_modify_new,stas_sp_current]=null_project_spatial_gpu(stas,mov,cell_params,matlab_cell_ids);
figure;
hist(mov_modify_new(:),50);
xlim([-200,200]);
violations = sum(abs(mov_modify_new(:))>maxClip+0.0001)
% distance = norm(mov_orig(:)-mov_modify_new(:))
z_k_half=2*mov_modify_new - mov;

x_k_1=z_k_half;
x_k_1(x_k_1>maxClip)=maxClip;
x_k_1(x_k_1<-maxClip)=-maxClip;

mov=mov + x_k_1 - mov_modify_new;

% Need to change post processing!!
togo = violations>0 % input('Continue Iterating?');
    end
   
    %mov_modify_new = gather(mov_modify_new);mov_orig=gather(mov_orig);
end

mov_new2 = mov_modify_new;


% output movie assignment
mov_new2 = mov_new2 * 0.5/127.5;
mov_orig2 = mov_orig2 * 0.5/127.5;

% get A

Filt_dim1=size(stas_sp_current{1},1);
Filt_dim2=size(stas_sp_current{1},2);
ncells=1;

A=zeros(ncells,Filt_dim1*Filt_dim2);
for icell=1:length(stas_sp_current)
    A(icell,:)=stas_sp_current{icell}(:)';
end

if(rank(A)<min(size(A)))
    display('Spatial STA matrix not well conditioned. Check selected cells.')
end

A = gpuArray(A);
[u,s,v]=svd(A,'econ');
Ainv=v*(s^-1)*u';

P_null = eye(size(Ainv,1),size(A,2))-(Ainv*A);

P_sta = (Ainv*A);

P_mat.P_null = P_null;
P_mat.P_sta = P_sta;