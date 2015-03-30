function [mov_orig,mov_modify_new]=null_space_movie2(datafile,cell_params,mov_params,solver)
%% Note: 
%cell_params.sta_spatial_method = 1 for just using 4th frame, 2 is for fitting spatial STA. 
% Use cell_params.use_fits=2 (clipping) if cell_params.sta_spatial_method = 1 and 
% use cell_params.use_fits=0 (no processing of STA) if
% cell_params.sta_spatial_method = 2;


%% Start parallel pool
if(0)
if ((exist('matlabpool')==2) && (matlabpool('size') == 0))
  try
    matlabpool open
  catch me
    warning('Failed to open parallel sessions using matlabpool:\n  %s\n',...
        me.message);
  end
end
end

if (exist('parpool')==2)
  try
    if (isempty(gcp('nocreate')))
      parpool
    end
  catch me
    warning('Failed to open parallel pool using parpool:\n  %s\n',...
        me.message);
  end
end


%% Load cells and STAs
global_vars
% global_vars2
%datafile = '2013-10-10-0/data000';
type_name= cell(1,1);
type_name{1}=cell_params.type_name_inp;

if(~strcmp(datafile,'load_from_cell_params'))
datarun=load_data(datafile)
datarun=load_sta(datarun)
datarun=load_params(datarun)

%get_cell_ids(datarun,type_name) % Vision IDs - used for loading fitted
%STAs!
if(strcmp(type_name{1},'userCellList'))
   idx=[1:length(datarun.cell_ids)];
   idx_list=[];
   for icell_check=1:length(cell_params.cell_list)
   idx_list=[idx_list;idx(datarun.cell_ids==cell_params.cell_list(icell_check))];
   end
   
   matlab_cell_ids=idx_list;
   clear idx_list
else
matlab_cell_ids=get_cell_indices(datarun,type_name);
end
stas=datarun.stas.stas(matlab_cell_ids);
else
stas=cell_params.stas;    
matlab_cell_ids=1;
end

% Load STAs

stas_new=cell(length(stas),1);
for icell=1:length(stas)
    st_temp=zeros(size(stas{1},2),size(stas{1},1),1,size(stas{1},4)); % DOUBT .. Could be a reason for things to fail!!!!!
    for itime=1:size(stas{1},4)
        st_temp(:,:,:,itime)=mean(stas{icell}(:,:,:,end-itime+1),3)'; % DOUBT .. Could be a reason for things to fail!!!!!
    end
    stas_new{icell}=st_temp;
end
stas=stas_new;
stas_orig=stas_new;

% Used in movie post process

[stas_clipped,totalMaskAccept,CellMasks]= clipSTAs(stas,cell_params);
mov_params.totalMaskAccept=totalMaskAccept;
cell_params.CellMasks=CellMasks;    
    
if(isfield(cell_params,'use_fits')==1)
    if(cell_params.use_fits==1)
    addpath(genpath('~/Nishal/matlab/code')); 
    addpath(genpath('~/Nishal/matlab/private/nishal/fwdfittingfunctions'));
    
    fit_info=cell(length(stas),1);
    display('Starting STA fitting');
    parfor issta=1:length(stas)
        issta
        fit_info{issta} = fit_sta(stas_new{issta});
    end
    
    
    full_fit=cell(1,1);
    cellSelected=zeros(length(stas),1);
    icnt=0;
   for issta=1:length(stas)
   if(fit_info{issta}~=[])
       icnt=icnt+1;
      full_fit{icnt} = sta_fit_function(fit_info{issta}.initial_params);
  cellSelected(issta)=1;
   else
   cellSelected(issta)=0;
   display(sprintf('Cell Removed %d',issta));
   end
    
   end
   stas=full_fit;
   mov_params.totalMaskAccept=ones(size(stas{1},1),size(stas{1},2));
   display('Using STA fits')
   
   % Run ClipSTA only for the mask.
      [stas_clipped,totalMaskAccept,CellMasks]= clipSTAs(stas,cell_params);
    mov_params.totalMaskAccept=totalMaskAccept;
    cell_params.CellMasks=CellMasks;    
    end
    
    if(cell_params.use_fits==2)
    [stas,totalMaskAccept,CellMasks]= clipSTAs(stas,cell_params);
    mov_params.totalMaskAccept=totalMaskAccept;
    stas_clipped=stas;
    cell_params.CellMasks=CellMasks;    
    end
    
n_cell=length(stas);
filt_len=size(stas{1},4);
var64=size(stas{1},1);
var32=size(stas{1},2);
filt_dim1=var64;
filt_dim2=var32;
%% Get /Generate Original movie
mov_params.var64=var64;
mov_params.var32=var32;

[mov,mov_params]=generate_movie(mov_params);
mov_orig=mov;
movie_time=mov_params.movie_time;
%% Solve ? 


% 'mov','mov_time' is the movie from previous part to be used here 


if(solver==1)
    
addpath('../lsqrSOL/');
damp=0;
atol=10^-6;
btol=10^-6;
conlim=1.0e+300; % Doubt!
itnlim=1000;
show=1;

b= -Ax(stas,mov,movie_time,n_cell);


tic;
[ x, istop, itn, r1norm, r2norm, Anorm, Acond, Arnorm, xnorm, var1 ]...
   = lsqrSOL( movie_time*n_cell, movie_time*filt_dim1*filt_dim2, @dual_AAtx2, b(:), damp, atol, btol, conlim, itnlim, show );
toc;

mov_modify_new=reshape(x,[filt_dim1,filt_dim2,movie_time])+mov;
res=norm(Ax(stas,mov_modify_new,movie_time,n_cell))
mov_modify_new = mov_modify_new*127.5/max(abs(mov_modify_new(:)));

end
if(solver == 2)
    
addpath('../craigSOL/');

atol=1.00e-06;%10^-8;
btol=1.00e-06;%10^-8;
conlim=1.00e+300; % Doubt!..   1.0e+12
itnlim=1000;
show=1;

b= -Ax(stas,mov,movie_time,n_cell);

tic;
[ x, y, istop, itn, rnorm, Anorm, Acond, xnorm ]...
   = craigSOL( movie_time*n_cell, movie_time*filt_dim1*filt_dim2, @dual_AAtx2, b(:), atol, btol, conlim, itnlim, show );
toc;

mov_modify_new=reshape(x,[filt_dim1,filt_dim2,movie_time])+mov;
res=norm(Ax(stas,mov_modify_new,movie_time,n_cell))
mov_modify_new = mov_modify_new*127.5/max(abs(mov_modify_new(:)));
    
end 

if(solver==3)
[~,mov_modify_new]=fourier_project(stas,mov);
end

if(solver==4)
[~,mov_modify_new]=null_project_spatial(stas,mov,cell_params,matlab_cell_ids);
end

if(solver==5) % 256x256
    [stas256x256,mov256x256]=preprocess256x256(stas,mov);
    [mov_orig256x256,mov_modify_new256x256]=fourier_project256x256(stas256x256,mov256x256);
    [mov_orig,mov_modify_new] =postprocess256x256(mov_orig256x256,mov_modify_new256x256,mov);
    
end

if(solver==6)
    for iter=1:5
        maxClip = (0.48/0.5)*127.5;
[~,mov_modify_new]=fourier_project(stas,mov);
violations = sum(abs(mov_modify_new(:))>maxClip)
mov_modify_new(mov_modify_new>maxClip)=maxClip;
mov_modify_new(mov_modify_new<-maxClip)=-maxClip;
mov=mov_modify_new;
norm(mov_orig(:)-mov_modify_new(:))
% Need to change post processing!!
    end
    
end


if(solver==7)
    for iter=1:5
            maxClip = (0.48/0.5)*127.5;
[~,mov_modify_new]=null_project_spatial(stas,mov,cell_params,matlab_cell_ids);
violations = sum(abs(mov_modify_new(:))>maxClip)
mov_modify_new(mov_modify_new>maxClip)=maxClip;
mov_modify_new(mov_modify_new<-maxClip)=-maxClip;
mov=mov_modify_new;
norm(mov_orig(:)-mov_modify_new(:))
% Need to change post processing!!
    end
    
end


if(solver==8) % Dykstra's Alternating Projections, Spatial
     togo=1;
    while togo==1
    maxClip = (0.48/0.5)*127.5;
[~,mov_modify_new]=null_project_spatial(stas,mov,cell_params,matlab_cell_ids);
figure;
hist(mov_modify_new(:),50);
xlim([-200,200]);
violations = sum(abs(mov_modify_new(:))>maxClip+0.0001)
distance = norm(mov_orig(:)-mov_modify_new(:))
z_k_half=2*mov_modify_new - mov;

x_k_1=z_k_half;
x_k_1(x_k_1>maxClip)=maxClip;
x_k_1(x_k_1<-maxClip)=-maxClip;

mov=mov + x_k_1 - mov_modify_new;

% Need to change post processing!!
togo = violations>0 % input('Continue Iterating?');
    end
    
end

if(solver==9) % Dykstra's Alternating Projections, 64x32, Spatio-temporal
    togo=1;
    while togo==1
       
        maxClip = (0.48/0.5)*127.5;

        [~,mov_modify_new]=fourier_project(stas,mov);

        figure;
hist(mov_modify_new(:),50);
xlim([-200,200]);
violations = sum(abs(mov_modify_new(:))>maxClip+0.0001)
distance = norm(mov_orig(:)-mov_modify_new(:))
z_k_half=2*mov_modify_new - mov;

x_k_1=z_k_half;
x_k_1(x_k_1>maxClip)=maxClip;
x_k_1(x_k_1<-maxClip)=-maxClip;

mov=mov + x_k_1 - mov_modify_new;

% Need to change post processing!!
togo = input('Continue Iterating?');
    end
end
    
    
if(solver==10) % Dykstra's Alternating Projections,128x128, Spatio-temporal
    togo=1;
    while togo==1
       
        maxClip = (0.48/0.5)*127.5;

     

         [stas128x128,mov128x128]=preprocess128x128(stas,mov);
         [mov_orig128x128,mov_modify_new128x128]=fourier_project128x128(stas128x128,mov128x128);
         [~,mov_modify_new] =postprocess128x128(mov_orig128x128,mov_modify_new128x128,mov);
  
    
        figure;
hist(mov_modify_new(:),50);
xlim([-200,200]);
violations = sum(abs(mov_modify_new(:))>maxClip+0.0001)
distance = norm(mov_orig(:)-mov_modify_new(:))
z_k_half=2*mov_modify_new - mov;

x_k_1=z_k_half;
x_k_1(x_k_1>maxClip)=maxClip;
x_k_1(x_k_1<-maxClip)=-maxClip;

mov=mov + x_k_1 - mov_modify_new;

% Need to change post processing!!
togo = input('Continue Iterating?');
    end
    
end
    
if(solver==11) % Dykstra's Alternating Projections , 64x64, Spatio-temporal
    togo=1;
    while togo==1
       
        maxClip = (0.48/0.5)*127.5;


        [stas64x64,mov64x64]=preprocess64x64(stas,mov);
        [mov_orig64x64,mov_modify_new64x64]=fourier_project64x64(stas64x64,mov64x64);
        [~,mov_modify_new] =postprocess64x64(mov_orig64x64,mov_modify_new64x64,mov);
    
        figure;
hist(mov_modify_new(:),50);
xlim([-200,200]);
violations = sum(abs(mov_modify_new(:))>maxClip+0.0001)
distance = norm(mov_orig(:)-mov_modify_new(:))
z_k_half=2*mov_modify_new - mov;

x_k_1=z_k_half;
x_k_1(x_k_1>maxClip)=maxClip;
x_k_1(x_k_1<-maxClip)=-maxClip;

mov=mov + x_k_1 - mov_modify_new;

% Need to change post processing!!
togo = input('Continue Iterating?');
    end
end

if(solver==12) % Dykstra's Alternating Projections,256x256, Spatio-temporal
    togo=1;
    while togo==1
       
        maxClip = (0.48/0.5)*127.5;

     
    [stas256x256,mov256x256]=preprocess256x256(stas,mov);
    [mov_orig256x256,mov_modify_new256x256]=fourier_project256x256(stas256x256,mov256x256);
    [mov_orig,mov_modify_new] =postprocess256x256(mov_orig256x256,mov_modify_new256x256,mov);
    
    
        figure;
hist(mov_modify_new(:),50);
xlim([-200,200]);
violations = sum(abs(mov_modify_new(:))>maxClip+0.0001)
distance = norm(mov_orig(:)-mov_modify_new(:))
z_k_half=2*mov_modify_new - mov;

x_k_1=z_k_half;
x_k_1(x_k_1>maxClip)=maxClip;
x_k_1(x_k_1<-maxClip)=-maxClip;

mov=mov + x_k_1 - mov_modify_new;

% Need to change post processing!!
togo = input('Continue Iterating?');
    end
    
end

if(solver==14) % Dykstra's Alternating Projections, Spatial, mean stuff for NSEM
     togo=1;
    while togo==1
    maxClip = (0.48/0.5)*127.5;
[~,mov_modify_new]=null_project_spatial(stas,mov,cell_params,matlab_cell_ids);
figure;
hist(mov_modify_new(:),50);
xlim([-200,200]);
mean(mov_modify_new(:))
violations = sum(abs(mov_modify_new(:)-mov_params.mean)>maxClip+0.0001)
distance = norm(mov_orig(:)-mov_modify_new(:))
z_k_half=2*mov_modify_new - mov;

x_k_1=z_k_half;
x_k_1(x_k_1>maxClip)=maxClip;
x_k_1(x_k_1<-maxClip)=-maxClip;

mov=mov + x_k_1 - mov_modify_new;

% Need to change post processing!!
togo = input('Continue Iterating?');
    end
    
end

%% see_movie
% see_movie2
%% Correct means ,etc ?? Movie correction left ? 
[mov_orig,mov_modify_new]=movie_post_process(mov_orig,mov_modify_new,mov_params);
see_movie2


%% Save movies and mask
save('/Volumes/Analysis/nishal/test-null.mat','mov_orig','mov_modify_new','mov_orig','CellMasks','totalMaskAccept');
end