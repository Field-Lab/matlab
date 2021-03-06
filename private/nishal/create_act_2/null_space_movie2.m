function [mov_orig,mov_modify_new]=null_space_movie2(datafile,cell_params,mov_params,solver)
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
%datafile = '2013-10-10-0/data000';
type_name= cell(1,1);
type_name{1}=cell_params.type_name_inp;

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
n_cell=length(stas);

% Load STAs

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
var64=64;
filt_dim1=var64;
filt_dim2=32;
%% Get /Generate Original movie
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

%% see_movie
see_movie2
%% Correct means ,etc ?? Movie correction left ? 
[mov_orig,mov_modify_new]=movie_post_process(mov_orig,mov_modify_new,mov_params);
see_movie2

mov_modify_new=mov_modify_new+mov_params.mean;
mov_orig=mov_orig+mov_params.mean;

end