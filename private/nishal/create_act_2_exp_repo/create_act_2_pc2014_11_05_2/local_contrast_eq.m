

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



%% Local constrast enhancement

% Get mov_log{1} ready

no_mov=10;
mov_ch=1;%30/no_mov;
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
mov_adapt_hist=zeros(32,64,movie_time);
for itime=1:movie_time
mov(:,:,itime)=(imresize(movie2(:,:,itime),[64,32])')/256; % Doubt!!% -127.5?
mov_adapt_hist(:,:,itime)=adapthisteq(mov(:,:,itime));
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
cell_resp_orig=Ax(stas,mov*256-128,movie_time,n_cell);
cell_resp_null=Ax(stas,mov_adapt_hist*256-128,movie_time,n_cell);


figure
subplot(3,1,1);
plot(cell_resp_orig);
ylim([-200,200]);
title('Original movie');



subplot(3,1,2);
plot(cell_resp_null);
ylim([-200,200])
title('Null movie response (bigger scale)');

subplot(3,1,3);
plot(cell_resp_null);
title('Null movie response (normal scale)');

%% Play movies
figure;

for itime=1:movie_time
    subplot(2,2,1);
    imagesc(mov(:,:,itime));
    axis image
    colorbar
    colormap(gray)
    title(sprintf('Time: %d',itime));
    
    subplot(2,2,2);
    imagesc(mov_adapt_hist(:,:,itime));
    axis image
    colorbar
    title(sprintf('Time: %d',itime));
    colormap(gray)
        
    subplot(2,2,3);
    imagesc(mov(:,:,itime)-mov_adapt_hist(:,:,itime));
    axis image
    colorbar
    title(sprintf('Time: %d',itime));
    colormap(gray)
    
    pause(1/240)

end

