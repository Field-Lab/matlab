%% 
% Code for loading STAs and stuff..
 matlabpool
%%
% param_path = '/Volumes/Analysis/2013-10-10-0/data000/data000.params';
% paramFile=edu.ucsc.neurobiology.vision.io.ParametersFile(param_path);
% % 
% %     'x0'                   'Double'     
% %     'y0'                   'Double'     
% %     'SigmaX'               'Double'     
% %     'SigmaY'               'Double'     
% %     'Theta'                'Double'
% % 
% %     'RedTimeCourse'        'DoubleArray'
% %     'GreenTimeCourse'      'DoubleArray'
% %     'BlueTimeCourse'       'DoubleArray'
% neuronID=121
%  paramFile.getDoubleCell(neuronID, 'x0')
% paramFile.getArrayCell(neuronID, 'RedTimeCourse')
% paramFile.getIDList()
% %paramFile.getNeuronsInClass('On Parasol') -> Not working

%%
global_vars
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




%% 
% make STA B/W - if not doing BW, still have to make it reverse .. as STA
% read is not reversed!!

stas_new=cell(length(stas),1);
for icell=1:length(stas)
    st_temp=zeros(size(stas{1},1),size(stas{1},2),1,size(stas{1},4));
    for itime=1:30
        st_temp(:,:,:,itime)=mean(stas{icell}(:,:,:,end-itime+1),3);
    end
    stas_new{icell}=st_temp;
end
stas=stas_new;
filt_len=size(stas{1},4);
%% play STA
figure
for itime=1:30
%imagesc(rgb2gray(stas{10}(:,:,:,end-itime+1))),colormap(gray)
imagesc(stas{10}(:,:,:,itime));
colormap(gray);
caxis([-0.03,0.1])
title(sprintf('Time %0.02f',itime))
axis image
pause(0.1)
end


%% Movie stuff
no_mov=10;
mov_ch=30/no_mov;
mov_log=cell(mov_ch,1);
% See movie


movie_time=120*no_mov;%size(movie,3);

for imov_chunk=1:mov_ch
    imov_chunk
movie_log=zeros(320,160,movie_time);
icnt=0;
for imov=(imov_chunk-1)*no_mov+1:(imov_chunk-1)*no_mov+no_mov
    imov
    icnt=icnt+1;
load(sprintf('/Volumes/Data/stimuli/movies/eye-movement/NS_brownian/nishal_test_chunks/movie_chunk_%d',imov));
movie_log(:,:,(icnt-1)*120+1:icnt*120)=movie;
end
movie2=movie_log;
clear movie_log

%mov = rand(filt_dim,filt_dim,movie_time);
%mov=movie/max(movie(:));
%mov=mov-0.5;%mean(mov(:));
mov=zeros(32,64,movie_time);
for itime=1:movie_time
mov(:,:,itime)=(imresize(movie2(:,:,itime),[64,32])'-128); % Doubt!!
end
% figure
% for itime=1:movie_time
%     subplot(2,1,1);
%     imagesc(movie(:,:,itime)');
%     caxis([0,255]);
%     axis image
%     subplot(2,1,2);
%     mov(:,:,itime)=(imresize(movie(:,:,itime),[64,32])'-128); % Doubt!!
%     % Better downsampling!?
%    imagesc(mov(:,:,itime))
%    caxis([-128,128]);
%    colormap gray
%    axis image
%    title(sprintf('Time: %d',itime));
%    %pause(1/120);
% end
mov_log{imov_chunk}=mov;

end

filt_dim1=32;
filt_dim2=64;


%%
for imov_chunk=1:mov_ch
grad_in_dual_need_for_speed(stas,mov_log{imov_chunk},movie_time,n_cell,filt_dim1,filt_dim2,filt_len,imov_chunk)
end
%%
mov_reshape=zeros(filt_dim1*filt_dim2*movie_time,1);
itime=0;
for i=1:filt_dim1*filt_dim2:length(mov_reshape)
    itime=itime+1;
    xx=mov(:,:,itime);
  mov_reshape(i:i+filt_dim1*filt_dim2-1)=  xx(:);
end

%%
% tic;
% filt_len=size(stas{1},4);
% 
% % making sparse toeplitz ?
% %filt_mat= sparse([],[],[],movie_time*n_cell,filt_dim1*filt_dim2*movie_time,n_cell*movie_time*6*320*160);
% %filt_mat= sparse([],[],[],movie_time*1,filt_dim1*filt_dim2*movie_time,n_cell*movie_time*6*320*160);
% filt_mat=cell(n_cell,1);
% inv_filt_mat=cell(n_cell,1);
% 
% for icell=1:n_cell
% k=stas{icell};
% nzero=1*movie_time*filt_len*filt_dim1*filt_dim2;
% x_log=zeros(nzero,1);
% y_log=zeros(nzero,1);
% v_log=zeros(nzero,1);
% cnt=0;
% 
% icell
% for itime=1:movie_time
%    
%     for iitime=1:min(itime,filt_len)
%     k_t=k(:,:,iitime);   
%     k_t=k_t(:)';
%     %filt_mat(itime+(icell-1)*(movie_time),itime*filt_dim1*filt_dim2-(iitime)*filt_dim1*filt_dim2+1:itime*filt_dim1*filt_dim2-(iitime-1)*filt_dim1*filt_dim2)=k_t(:)';
%     for ispace=1:filt_dim1*filt_dim2
%         cnt=cnt+1;
%     x_log(cnt)=itime+(1-1)*(movie_time);
%     y_log(cnt)=itime*filt_dim1*filt_dim2-(iitime)*filt_dim1*filt_dim2+ispace;
%     v_log(cnt)=k_t(ispace);
%         
%     end
%     
%     end
% end
% x_log=x_log(1:cnt);
% y_log=y_log(1:cnt);
% v_log=v_log(1:cnt);
% 
% filt_mat{icell}=sparse(x_log,y_log,v_log,movie_time*1,filt_dim1*filt_dim2*movie_time,1*movie_time*6*320*160);
% inv_filt_mat{icell}=(filt_mat{icell}*filt_mat{icell}')\filt_mat{icell};
% end
% toc;
%%

% Parallel version

tic;
filt_len=size(stas{1},4);

% making sparse toeplitz ?
%filt_mat= sparse([],[],[],movie_time*n_cell,filt_dim1*filt_dim2*movie_time,n_cell*movie_time*6*320*160);
%filt_mat= sparse([],[],[],movie_time*1,filt_dim1*filt_dim2*movie_time,n_cell*movie_time*6*320*160);
filt_mat=cell(n_cell,1);
inv_filt_mat=cell(n_cell,1);

parfor icell=1:n_cell % Does not use all cores, but faster than before!
    icell
[filt_mat{icell},inv_filt_mat{icell}]=make_filter_matrix(stas{icell},movie_time,filt_len,filt_dim1,filt_dim2);
end
toc;
