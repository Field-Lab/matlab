WN_datafile = '/Volumes/Analysis/2016-02-17-1/d28_51-norefit/data028/data028';
datarun=load_data(WN_datafile)
datarun=load_params(datarun);
datarun=load_sta(datarun);

vision_id=1801;
iidx=1:length(datarun.cell_ids);
mcellid = iidx(datarun.cell_ids==vision_id);
h=plot_rf_summaries(datarun, vision_id);hold on;

figure;
imagesc(sum(datarun.stas.stas{mcellid}(:,:,:,26),3));axis image;colormap gray;hold on;
XD = h.Children.XData;
YD = h.Children.YData;
plot(XD,YD,'Linewidth',2);
ylim([19,26]);
xlim([49,56]);
set(gca,'xTick',[]);
set(gca,'yTick',[]);

%%  Compute null movie for ONE cell

startup_bertha

startup_rooster

datafile='/Volumes/Analysis/2016-02-17-1/d28_51-norefit/data028/data028';%'nishal/2014-08-20-2/data001/data001';

% type_name_inp = 'userCellList' for a list of cells.

no_images_per_movie=10;
start_image=10;
tag='test-tag';
destination_mat=['/Volumes/Lab/Users/bhaishahster/',tag];
if(~exist(destination_mat,'dir'))
mkdir(destination_mat);
end

dest_raw='/Volumes/Data/stimuli/movies/null_space/';

movies=cell(20,1);



%% low contrast WN, spatial null
cell_params=struct();
cell_params.type_name_inp='userCellList';
cell_params.cell_list=vision_id;
cell_params.STAlen=14;
cell_params.sta_spatial=sprintf('%s/stas_spatial.mat',destination_mat);
cell_params.use_fits=2; % 2, 0,0,2
cell_params.sta_spatial_method=4;%1,2 ,3,4
% cell_params.sta_spatial_method = 1 for just using 4th frame, 2 is for fitting spatial STA. 
% Use cell_params.use_fits=2 (clipping) if cell_params.sta_spatial_method = 1 and 
% use cell_params.use_fits=0 (no processing of STA) if
% cell_params.sta_spatial_method = 2;
% STA spatial null Method 3 = low rank, 4 = average waveform and use it ..  

mov_params=struct();
mov_params.mov_type='bw-precomputed';
mov_params.movie_time=120*10/2;
mov_params.mean=0.5*255;
mov_params.deviation=0.24*255;

mov_params.mdf_file = '/Volumes/Analysis/stimuli/white-noise-xml/BW-8-2-0.48-11111.xml';
mov_params.stixel=8;

mov_params.interval = 2; % Not important to have this parameter. Default is 1. When we want repeated frames, just set this interval (This just controls the blank Frames), and select the movie_time appropriately.

% Post process. Default is stretch. If using default, need to give only mov_params.scaling_loss parameter.
mov_params.post_process_method = 'scale'; % or, 'stretch'
mov_params.scale = 0.48/0.48;
%mov_params.scaling_loss=0.05; % a number in [0,1], fraction of values that is changed by scaling.

solver=8; % Solver 4 used for spatial nulling, 7 for iterated spatial nulling
[mov_orignial,mov_modify_new]=null_space_movie2(datafile,cell_params,mov_params,solver);


movies{5}=mov_orignial;
movies{6}=mov_modify_new;
mov_idx=5;
write_movie_idx(destination_mat,movies{mov_idx},mov_idx,mov_params.stixel);
mov_idx=6;
write_movie_idx(destination_mat,movies{mov_idx},mov_idx,mov_params.stixel);
%%
% make movies
interval=4;
condMov=cell(3,1);
rawMovFrames=1200/(4);
icnt=0;
% make pixel histogram
for imov=[5,6]
    %[stim,height,width,header_size] = get_raw_movie(sprintf('/Volumes/Lab/Users/bhaishahster/test-tag/%d.rawMovie',imov),rawMovFrames,1);
    [stim,height,width,header_size] = get_raw_movie(sprintf('/Volumes/Data/2016-02-17-1/Visual/null/pc2016_02_17_1_data028/%d.rawMovie',imov),rawMovFrames,1);
  
    subtract_movies{3}=mean(stim,1);
    subtract_movies{3}=mean(stim,1)*0+127.5;
    movie=stim-repmat(subtract_movies{3},[rawMovFrames,1,1]);
    movie=movie/255;
    
    icnt=icnt+1;
    qq=permute(movie,[2,3,1]);
    ifcnt = 0;
    condMov{icnt}=zeros(size(qq,1),size(qq,2),size(qq,3)*interval);
    for iframe=1:size(qq,3)
        for irepeat=1:interval
            ifcnt=ifcnt+1;
            condMov{icnt}(:,:,ifcnt)=qq(:,:,iframe)+0.5; % cond mov is between 0 and 1 now!
        end
        
    end
    
    condMov{icnt}=permute(condMov{icnt},[2,1,3]);
end

iframe=457;
figure;
imagesc(condMov{1}(:,:,iframe));axis image;colormap gray;hold on;
plot(XD,YD,'Linewidth',2);
ylim([19,26]);
xlim([49,56]);
set(gca,'xTick',[]);
set(gca,'yTick',[]);
caxis([0,1]);

figure;
imagesc(condMov{2}(:,:,iframe));axis image;colormap gray;hold on;
plot(XD,YD,'Linewidth',2);
ylim([19,26]);
xlim([49,56]);
set(gca,'xTick',[]);
set(gca,'yTick',[]);
caxis([0,1]);


