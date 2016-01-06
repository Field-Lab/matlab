addpath(genpath('../null_analyse/'));
addpath(genpath('../null_analyse/analyse_functions'));
%startup_null_analyse_tenessee
startup_null_analyse_bertha

%%

% Condition strings
nConditions=8;
condDuration=1200/4;
cond_str=cell(nConditions,1);
cond_str{1}='Original';
cond_str{2}='Null ';

interestingConditions=[1,2,3,4,5,6,7,8];

dataRuns_OFF_additivity = [23,24,26,27,29,31,34,35];
dataRuns_ON_additivity = [23,25,26,28,30,32,34,36];
movies_OFF_addivitiy =[1,2,3,4,14,18,5,6];
movies_ON_additivity = [1,8,3,10,16,20,5,12];

%% Load Movies
%  
rawMovFrames=1200/(4);
figure;
icnt=0;
% make pixel histogram
for imov=movies_ON_additivity
    [stim,height,width,header_size] = get_raw_movie(sprintf('/Volumes/Data/2015-11-09-8/Visual/2015-11-09-8-from007/%d.rawMovie',imov),rawMovFrames,1);
    subtract_movies{3}=mean(stim,1);
    subtract_movies{3}=mean(stim,1)*0+127.5;
    movie=stim-repmat(subtract_movies{3},[rawMovFrames,1,1]);
    movie=movie/255;
    
    icnt=icnt+1;
    
    subplot(4,2,icnt);
    qq=movie;
    hist(qq(:),20)
    title(sprintf('Movie pixel histogram %d',imov));
end

% % make movies
% interval=4;
% condMov=cell(nConditions,1);
% rawMovFrames=1200/(4);
% icnt=0;
% % make pixel histogram
% for imov=[1,2,4,10,12,5,6,8]
%     [stim,height,width,header_size] = get_raw_movie(sprintf('/Volumes/Data/2015-09-23-0/Visual/2015-09-23-0_data017/%d.rawMovie',imov),rawMovFrames,1);
%     subtract_movies{3}=mean(stim,1);
%     subtract_movies{3}=mean(stim,1)*0+127.5;
%     movie=stim-repmat(subtract_movies{3},[rawMovFrames,1,1]);
%     movie=movie/255;
%     
%     icnt=icnt+1;
%     qq=permute(movie,[2,3,1]);
%     ifcnt = 0;
%     condMov{icnt}=zeros(size(qq,1),size(qq,2),size(qq,3)*interval);
%     for iframe=1:size(qq,3)
%         for irepeat=1:interval
%             ifcnt=ifcnt+1;
%             condMov{icnt}(:,:,ifcnt)=qq(:,:,iframe)+0.5; % cond mov is between 0 and 1 now!
%         end
%         
%     end
%     
% end
% 
% make contrast map
rawMovFrames=1200/(4);
figure;
icnt=0;
cMap = cell(8,1);
h=figure('Color','w');
for imov=movies_ON_additivity
    [stim,height,width,header_size] = get_raw_movie(sprintf('/Volumes/Data/2015-11-09-8/Visual/2015-11-09-8-from007/%d.rawMovie',imov),rawMovFrames,1);
    subtract_movies{3}=mean(stim,1);
    subtract_movies{3}=mean(stim,1)*0+127.5;
    movie=stim-repmat(subtract_movies{3},[rawMovFrames,1,1]);
    movie=movie/255;
    
    icnt=icnt+1;
    
    subplot(3,2,icnt);
    qq=movie;
    
    cMap{icnt}=contrastMap(qq);
    if(icnt ==1 ||icnt==3)
    imagesc(cMap{icnt}');
    title(sprintf('cMap: %d',icnt));
    colormap gray
    caxis([0,1]);
    end
    
    if(icnt==2)
    imagesc(( cMap{icnt}-cMap{1})'./cMap{1}');
    title(sprintf('cMap: %d, compared to 1',icnt));
    caxis([-1,1]);
    colormap gray
    end
    
    if(icnt==4 || icnt==5 ||icnt==6)
    imagesc(( cMap{icnt}-cMap{3} )'./cMap{3}');
        title(sprintf('cMap: %d, compared to 3',icnt));
        caxis([-1,1]);
        colormap gray
    end
   % caxis([3,6]);
    colorbar
    axis image
 
end

% s=hgexport('readstyle','cMap');
% hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_09_23_0/data000/cMap.eps'),s);



%% additivity experiment -  OFF

dataRuns =dataRuns_ON_additivity;

WN_datafile = '/Volumes/Analysis/2015-12-18-2/d22_37-norefit/data015/data015';


datarun=load_data(WN_datafile)
datarun=load_params(datarun)

cellTypeId = 1
InterestingCell_vis_id=datarun.cell_types{cellTypeId}.cell_ids; 
cellTypeUsed=cellTypeId*ones(length(InterestingCell_vis_id),1);

condDuration=10;
nConditions=1;

cols='rkrkrkrkrkrkkrkrkrkr';
spkCondColl=cell(length(dataRuns),1);

for ref_cell_number=1:length(InterestingCell_vis_id); %11
    close all
    cellID= InterestingCell_vis_id(ref_cell_number)
        
    for idata=1:length(dataRuns)
    Null_datafile = sprintf('/Volumes/Analysis/2015-12-18-2/d22_37-norefit/data0%02d',dataRuns(idata));
    neuronPath = [Null_datafile,sprintf('/data0%02d.neurons',dataRuns(idata))];
    [spkColl,spkCondColl{idata},h]=plot_raster_script_pc2015_12_18_2_light_photons(cellID,nConditions,condDuration,cond_str,neuronPath);
    end
    
    
    h=figure;
    for irun = 1:length(dataRuns)
    plot(spkCondColl{irun}.xPoints/20000,spkCondColl{irun}.yPoints - (irun-1)*30,cols(irun));
    hold on;
    end
    set(gca,'yTick',[]);
    ylim([-(length(dataRuns)-1)*30,30]);
    InterestingCell_vis_id(ref_cell_number)
    pause(1/120);
     if(~isdir(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_12_18_2/d_add/CellType_%s',datarun.cell_types{cellTypeId}.name)))
         mkdir(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_12_18_2/d_add/CellType_%s',datarun.cell_types{cellTypeId}.name));
     end
    s=hgexport('readstyle','raster');
    hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_12_18_2/d_add/CellType_%s/CellID_%d.eps',datarun.cell_types{cellTypeId}.name,InterestingCell_vis_id(ref_cell_number)),s);
end

%% Long null
%% long null experiment -  ON

dataRuns =28%
imov = 20%20,18;

% Load movie
interval=4;
rawMovFrames=216000/(interval);
figure;
icnt=0;
[stim,height,width,header_size] = get_raw_movie(sprintf('/Volumes/Data/2015-11-09-8/Visual/2015-11-09-8-from007/%d.rawMovie',imov),rawMovFrames,1);
subtract_movies{3}=mean(stim,1);
subtract_movies{3}=mean(stim,1)*0+127.5;
movie=stim-repmat(subtract_movies{3},[rawMovFrames,1,1]);
movie=movie/255;

movie = permute(movie,[2,3,1]);
movie=repelem(movie,1,1,interval);


% Load spikes
WN_datafile = '/Volumes/Analysis/2015-11-09-8/d15_61-norefit/data007/data007';

datarun=load_data(WN_datafile)
datarun=load_params(datarun)
datarun = load_sta(datarun);

cellTypeId = 1;
InterestingCell_vis_id=[2926,2971,2854,5642];%datarun.cell_types{cellTypeId}.cell_ids; 
cellTypeUsed=cellTypeId*ones(length(InterestingCell_vis_id),1);

condDuration=1800;
nConditions=1;

cols='rkrkrkrkrkrkkrkrkrkr';
spkCondColl=cell(length(dataRuns),1);
close all
for ref_cell_number = 1:length(InterestingCell_vis_id)
 cellID=InterestingCell_vis_id(ref_cell_number);
        
    idata=1;
    Null_datafile = sprintf('/Volumes/Analysis/2015-11-09-8/d15_61-norefit/data0%02d',dataRuns(idata));
    neuronPath = [Null_datafile,sprintf('/data0%02d.neurons',dataRuns(idata))];
   [spks,spksGen_hr]=get_spikes_pc2015_11_09_1(neuronPath,cellID,condDuration);
    
    staDepth=30;
    STA = calc_sta_mov(movie,spks',staDepth);
%     figure;
%     for itime=1:staDepth%1:staDepth
%     imagesc(STA(:,:,itime)');
%     title(sprintf('%d',itime));
%     colormap gray
%     axis image
%     caxis([min(STA(:)), max(STA(:))]);
%     pause;
%    end
    
   %end

 iidx = 1:length(datarun.cell_ids); 
 matlabID = iidx(cellID == datarun.cell_ids);
 WNsta = datarun.stas.stas{matlabID};
% figure;
% imagesc(squeeze(mean(WNsta(:,:,:,26),3)));axis image;colormap gray

icnt=1;
scrsz = get(groot,'ScreenSize');
h=figure('Position', [1 scrsz(4) scrsz(3) scrsz(4)]);
for itime = [4,6,8,10,12,14,16,18,20,22]
    
    subplot(2,10,icnt);
    imagesc(repelem(STA(:,:,itime),20,20));
    title(sprintf('%d',itime));
    colormap gray
    axis image
    caxis([min(STA(:)), max(STA(:))]);
    icnt=icnt+1;
    title(sprintf('%d',itime),'FontSize',10);
    set(gca,'xTick',[]); set(gca,'yTick',[]);
end
for itime = [4,6,8,10,12,14,16,18,20,22]
    subplot(2,10,icnt);
    imagesc(repelem(squeeze(mean(WNsta(:,:,:,end-itime+1),3))',20,20));
    title(sprintf('%d',itime));
    colormap gray
    axis image
    caxis([min(WNsta(:)), max(WNsta(:))]);
    icnt=icnt+1;
    title(sprintf('%d',itime),'FontSize',10);
    set(gca,'xTick',[]); set(gca,'yTick',[]);
     
end
suptitle(sprintf('Cell ID: %d',cellID));

   pause(1/120);
     if(~isdir(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_11_09_8/d_long_STA/CellType_%s',datarun.cell_types{cellTypeId}.name)))
         mkdir(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_11_09_8/d_long_STA/CellType_%s',datarun.cell_types{cellTypeId}.name));
     end
    s=hgexport('readstyle','raster');
    hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_11_09_8/d_long_STA/CellType_%s/CellID_%d.eps',datarun.cell_types{cellTypeId}.name,cellID),s);
    end
%% Analyse a cell - fit ASM, difference of rasters, non-linearitiy ..

%% 
WN_datafile = '/Volumes/Analysis/2015-11-09-8/d15_61-norefit/data007/data007';
datarun=load_data(WN_datafile)
datarun=load_params(datarun)


cellTypeId = 1%OFF:2 .. 
cellID_select=[2926,2971,2854,5642]% datarun.cell_types{cellTypeId}.cell_ids; % 51 ? [6741,2176,6961]


 %% fit ASM
WN_datafile = '/Volumes/Analysis/2015-11-09-8/d15_61-norefit/data007/data007';

movie_xml = 'RGB-8-2-0.48-11111';
stim_length=900;% 
cellID_list = [cellID_select];%1006,1411,2086,2341,2686,4096,4276,4831,5371,5566,5596,5866,7156,7216,7381,7490];
nSU_list= [1,2,3,4,5,6,7,8,9,10];
destination= 'pc2015_11_09_8_analysis_fits/SUs_data007'
sta_depth=30; % CAREFUL TO CHANGE IT TO LONG .. FOR SLOWER CELLS

[~] = get_gmlm_sta2(WN_datafile,movie_xml,stim_length,cellID_list,nSU_list,destination,sta_depth)

%% Learn from low contrast WN in data 029

 %% fit ASM -mask from WN movie, but fit to custom input and spikes
 % load long null movie
 imov = 20;

% Load movie
interval=4;
rawMovFrames=216000/(interval);
figure;
icnt=0;
[stim,height,width,header_size] = get_raw_movie(sprintf('/Volumes/Data/2015-11-09-8/Visual/2015-11-09-8-from007/%d.rawMovie',imov),rawMovFrames,1);
subtract_movies{3}=mean(stim,1);
subtract_movies{3}=mean(stim,1)*0+127.5;
movie=stim-repmat(subtract_movies{3},[rawMovFrames,1,1]);
movie=movie/255;

movie = permute(movie,[2,3,1]);
movie=repelem(movie,1,1,interval);


 % WN details to get the significant stixels box   
WN_datafile = '/Volumes/Analysis/2015-11-09-8/d15_61-norefit/data007/data007';
movie_xml = 'RGB-8-2-0.48-11111';
stim_length=1800;% 
% load cell IDs for ON type 

% Load spikes
datarun=load_data(WN_datafile)
datarun=load_params(datarun)
datarun = load_sta(datarun);
cellTypeId = 1;
InterestingCell_vis_id=[2854,5642]%datarun.cell_types{cellTypeId}.cell_ids; 
idata=1;

for cellID_list = InterestingCell_vis_id;
    
 % spikes for null movie
condDuration=1800;dataRuns = 28;% 25 - ON , 27 - OFF
Null_datafile = sprintf('/Volumes/Analysis/2015-11-09-8/d15_61-norefit/data0%02d',dataRuns(idata));
neuronPath = [Null_datafile,sprintf('/data0%02d.neurons',dataRuns(idata))];
[spks,spksGen_hr]=get_spikes_pc2015_11_09_1(neuronPath,cellID_list,condDuration);
   
nSU_list= [1,2,3,4,5,6,7,8,9,10];
destination= 'pc2015_11_09_8_analysis_fits/SUs_data028_exp'
sta_depth=30; % CAREFUL TO CHANGE IT TO LONG .. FOR SLOWER CELLS

[~] = get_gmlm_sta4(WN_datafile,movie_xml,stim_length,cellID_list,nSU_list,destination,sta_depth,movie,spks)
end

