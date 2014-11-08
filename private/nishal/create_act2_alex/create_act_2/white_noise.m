function [mov]=white_noise(datafile,type_name_inp,destination_mat,movie_white_len,earliest_peak_time)


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
%%
global_vars
%datafile = '2013-10-10-0/data000';
type_name= cell(1,1);
type_name{1}=type_name_inp;

datarun=load_data(datafile)
datarun=load_sta(datarun)
datarun=load_params(datarun)

%get_cell_ids(datarun,type_name) % Vision IDs - used for loading fitted
%STAs!

matlab_cell_ids=get_cell_indices(datarun,type_name);
stas=datarun.stas.stas(matlab_cell_ids);
n_cell=length(stas);


%% Load STAS

stas_new=cell(length(stas),1);
for icell=1:length(stas)
    st_temp=zeros(size(stas{1},2),size(stas{1},1),1,size(stas{1},4)); % DOUBT .. Could be a reason for things to fail!!!!!
    for itime=1:30
        st_temp(:,:,:,itime)=mean(stas{icell}(:,:,:,end-itime+1),3)'; % DOUBT .. Could be a reason for things to fail!!!!!
    end
    stas_new{icell}=st_temp;
end
stas=stas_new;
filt_len=size(stas{1},4);



%% Make whitish noise movie

var64=64;
filt_dim1=var64;
filt_dim2=32;
mov=zeros(filt_dim1,filt_dim2,movie_white_len);
init_mov=zeros(filt_dim1,filt_dim2,movie_white_len);
A_mat = zeros(n_cell,filt_dim1*filt_dim2);
for icell=1:n_cell
x=stas{icell}(:,:,1,earliest_peak_time);
A_mat(icell,:)=x(:)';
end
clear x
tic;
A_right_inv = (A_mat'*((A_mat*A_mat')\eye(n_cell,n_cell)));
toc;

tic;
for itime=filt_len:movie_white_len
   
    
    tot_resp=zeros(n_cell,1);
   
    mov_chunk=mov(:,:,itime-filt_len+earliest_peak_time:itime-1);
    for icell=1:n_cell
       tot_resp(icell)= sum(sum(sum(squeeze(stas{icell}(:,:,1,end:-1:earliest_peak_time+1)).*mov_chunk)));
    end
    
    init_frame= randn(filt_dim1,filt_dim2);
    b = - A_mat*init_frame(:) - tot_resp;
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

%% Play mov
figure;
for itime=1:120:movie_white_len
    subplot(2,1,1);
    imagesc(init_mov(30:35,10:15,itime));
    colormap gray
    colorbar
    
    subplot(2,1,2);
    imagesc(mov(30:35,10:15,itime));
    colormap gray
    colorbar
    %caxis([-1,2]);
    %title(sprintf('Time: %d',itime));
    pause(1);
end
end
