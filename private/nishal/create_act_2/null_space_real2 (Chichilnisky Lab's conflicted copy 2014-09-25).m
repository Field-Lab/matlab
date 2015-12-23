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

%% play STA
figure
for itime=1:30
%imagesc(rgb2gray(stas{10}(:,:,:,end-itime+1))),colormap(gray)
imagesc(stas{1}(:,:,:,itime));
colormap(gray);
caxis([-0.03,0.1])
title(sprintf('Time %0.02f',itime))
pause(1)

end


%% Movie stuff

for imov=6:6
load(sprintf('/Volumes/Data/stimuli/movies/eye-movement/NS_brownian/nishal_test_chunks/movie_chunk_%d',imov));
end


% See movie
movie_time=120;%size(movie,3);
%mov = rand(filt_dim,filt_dim,movie_time);
%mov=movie/max(movie(:));
%mov=mov-0.5;%mean(mov(:));
mov=zeros(32,64,movie_time);
figure
for itime=1:movie_time
    subplot(2,1,1);
    imagesc(movie(:,:,itime)');
    caxis([0,255]);
    axis image
    subplot(2,1,2);
    mov(:,:,itime)=imresize(movie(:,:,itime),[64,32])'-128; 
    % Better downsampling!?
   imagesc(mov(:,:,itime))
   caxis([-128,128]);
   colormap gray
   axis image
   title(sprintf('Time: %d',itime));
   pause(1/120);
end

filt_dim1=32;
filt_dim2=64;
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
