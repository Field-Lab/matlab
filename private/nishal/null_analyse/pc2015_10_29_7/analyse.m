addpath(genpath('../null_analyse/'));
addpath(genpath('../null_analyse/analyse_functions'));
%startup_null_analyse_tenessee
startup_null_analyse_bertha

%%
% Condition strings
nConditions=6;
condDuration=1200/4;
cond_str=cell(nConditions,1);
cond_str{1}='Original';
cond_str{2}='Null ';

interestingConditions=[1,2,3,4,5,6,7];

dataRuns_OFF_additivity = [14,15,17,18,20,22,24];
dataRuns_ON_additivity = [14,16,17,19,21,23,24];
movies_OFF_addivitiy =[1,2,5,6,10,14,11];
movies_ON_additivity = [1,4,5,8,12,16,11];
%% Load Movies
%  
rawMovFrames=1200/(4);
figure;
icnt=0;
% make pixel histogram
for imov=movies_ON_additivity
    [stim,height,width,header_size] = get_raw_movie(sprintf('/Volumes/Data/2015-10-29-7/Visual/%d.rawMovie',imov),rawMovFrames,1);
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
%     [stim,height,width,header_size] = get_raw_movie(sprintf('/Volumes/Data/2015-11-09-1/Visual/pc2015_11_09_1_fromdata009/%d.rawMovie',imov),rawMovFrames,1);
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
%% contrast map 
% set up stuff for reading real responses

WN_datafile ='/Volumes/Analysis/2015-10-29-7/streamed/data012/data012'

datarun=load_data(WN_datafile)
datarun=load_params(datarun)
cellTypeId = 1;
InterestingCell_vis_id=datarun.cell_types{cellTypeId}.cell_ids; 


% make contrast map
rawMovFrames=1200/(4);
figure;
icnt=0;
cMap = cell(8,1);
h=figure('Color','w');
for imov=movies_ON_additivity
    [stim,height,width,header_size] = get_raw_movie(sprintf('/Volumes/Data/2015-10-29-7/Visual/%d.rawMovie',imov),rawMovFrames,1);
    subtract_movies{3}=mean(stim,1);
    subtract_movies{3}=mean(stim,1)*0+127.5;
    movie=stim-repmat(subtract_movies{3},[rawMovFrames,1,1]);
    movie=movie/255;
    
    icnt=icnt+1;
    
    subplot(3,2,icnt);
    qq=movie;
    
    cMap{icnt}=contrastMap(qq);
    if(icnt ==1 ||icnt==3)
        imagesc(repelem(cMap{icnt}(:,end:-1:1)',20,20));
        title(sprintf('cMap: %d',icnt));
        colormap gray;axis image
        caxis([0,0.5]);set(gca,'xTick',[]);set(gca,'yTick',[]);
        plot_rf_fit_nishal(datarun, InterestingCell_vis_id,'magnify',20,'fill_color',[0,1,1],'alpha',0.3,'fill',false,'edge',true,'labels',false,'edge_color',[1,0,0]);
    end
    
    if(icnt==2)
  
        h=imagesc(repelem(( cMap{icnt}(:,end:-1:1)-cMap{1}(:,end:-1:1))'./cMap{1}(:,end:-1:1)',20,20));
        title(sprintf('cMap: %d, compared to 1',icnt));
        caxis([-1,1]);set(gca,'xTick',[]);set(gca,'yTick',[]);
        colormap gray;axis image
        plot_rf_fit_nishal(datarun, InterestingCell_vis_id,'magnify',20,'fill_color',[0,1,1],'alpha',0.3,'fill',false,'edge',true,'labels',false,'edge_color',[1,0,0]);
    end
    
    if(icnt==4 || icnt==5 ||icnt==6)
        imagesc(repelem(( cMap{icnt}(:,end:-1:1)-cMap{3}(:,end:-1:1) )'./cMap{3}(:,end:-1:1)',20,20));
        title(sprintf('cMap: %d, compared to 3',icnt));
        caxis([-1,1]);axis image
        colormap gray;set(gca,'xTick',[]);set(gca,'yTick',[]);
        plot_rf_fit_nishal(datarun, InterestingCell_vis_id,'magnify',20,'fill_color',[0,1,1],'alpha',0.3,'fill',false,'edge',true,'labels',false,'edge_color',[1,0,0]);
    end
    % caxis([3,6]);
    colorbar
    axis image
    
end

% s=hgexport('readstyle','cMap');
% hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_09_23_0/data000/cMap.eps'),s);

mask = (cMap{2}-cMap{1})~=0;
imagesc(mask);


%% additivity experiment 

dataRuns =dataRuns_OFF_additivity;

WN_datafile = '/Volumes/Analysis/2015-10-29-7/d12_41-norefit/data012/data012';


datarun=load_data(WN_datafile)
datarun=load_params(datarun)

cellTypeId = 2
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
    Null_datafile = sprintf('/Volumes/Analysis/2015-10-29-7/d12_41-norefit/data0%02d',dataRuns(idata));
    neuronPath = [Null_datafile,sprintf('/data0%02d.neurons',dataRuns(idata))];
    [spkColl,spkCondColl{idata},h]=plot_raster_script_pc2015_09_23_0_light(cellID,nConditions,condDuration,cond_str,neuronPath);
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
     if(~isdir(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_10_29_7/d_add/CellType_%s',datarun.cell_types{cellTypeId}.name)))
         mkdir(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_10_29_7/d_add/CellType_%s',datarun.cell_types{cellTypeId}.name));
     end
    s=hgexport('readstyle','raster');
    hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_10_29_7/d_add/CellType_%s/CellID_%d.eps',datarun.cell_types{cellTypeId}.name,InterestingCell_vis_id(ref_cell_number)),s);
end

%% Long null
%% long null experiment -  ON

dataRuns =25%27,25;
imov = 18%20,18;

% Load movie
interval=4;
rawMovFrames=216000/(interval);
figure;
icnt=0;
[stim,height,width,header_size] = get_raw_movie(sprintf('/Volumes/Data/2015-11-09-1/Visual/pc2015_11_09_1_fromdata009/%d.rawMovie',imov),rawMovFrames,1);
subtract_movies{3}=mean(stim,1);
subtract_movies{3}=mean(stim,1)*0+127.5;
movie=stim-repmat(subtract_movies{3},[rawMovFrames,1,1]);
movie=movie/255;

movie = permute(movie,[2,3,1]);
movie=repelem(movie,1,1,interval);


% Load spikes
WN_datafile = '/Volumes/Analysis/2015-11-09-1/d14_36-norefit/data009/data009';

datarun=load_data(WN_datafile)
datarun=load_params(datarun)
datarun = load_sta(datarun);

cellTypeId = 1;
InterestingCell_vis_id=datarun.cell_types{cellTypeId}.cell_ids; 
cellTypeUsed=cellTypeId*ones(length(InterestingCell_vis_id),1);

condDuration=1800;
nConditions=1;

cols='rkrkrkrkrkrkkrkrkrkr';
spkCondColl=cell(length(dataRuns),1);
close all
for ref_cell_number = 1:length(InterestingCell_vis_id)
 cellID=InterestingCell_vis_id(ref_cell_number);
        
    idata=1;
    Null_datafile = sprintf('/Volumes/Analysis/2015-11-09-1/d14_36-norefit/data0%02d',dataRuns(idata));
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
     if(~isdir(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_11_09_1/d_long_STA/CellType_%s',datarun.cell_types{cellTypeId}.name)))
         mkdir(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_11_09_1/d_long_STA/CellType_%s',datarun.cell_types{cellTypeId}.name));
     end
    s=hgexport('readstyle','raster');
    hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_11_09_1/d_long_STA/CellType_%s/CellID_%d.eps',datarun.cell_types{cellTypeId}.name,cellID),s);
    end
%% Analyse a cell - fit ASM, difference of rasters, non-linearitiy ..

%% 
WN_datafile = '/Volumes/Analysis/2015-11-09-1/d14_36-norefit/data009/data009';
datarun=load_data(WN_datafile)
datarun=load_params(datarun)


cellTypeId = 1%OFF:2 .. 
cellID_select= datarun.cell_types{cellTypeId}.cell_ids; % 51 ? [6741,2176,6961]


 %% fit ASM
WN_datafile = '/Volumes/Analysis/2015-11-09-1/d14_36-norefit/data009/data009';

movie_xml = 'RGB-8-2-0.48-11111';
stim_length=1800;% 
cellID_list = 5341%[cellID_select];%1006,1411,2086,2341,2686,4096,4276,4831,5371,5566,5596,5866,7156,7216,7381,7490];
nSU_list= [1,2,3,4,5,6,7,8,9,10];
destination= 'pc2015_11_09_1_analysis_fits/SUs_data009'
sta_depth=30; % CAREFUL TO CHANGE IT TO LONG .. FOR SLOWER CELLS

[~] = get_gmlm_sta2(WN_datafile,movie_xml,stim_length,cellID_list,nSU_list,destination,sta_depth)

 %% fit ASM -mask from WN movie, but fit to custom input and spikes
 % load long null movie
 imov = 18%20 - OFF,18 - ON;

% Load movie
interval=4;
rawMovFrames=216000/(interval);
figure;
icnt=0;
[stim,height,width,header_size] = get_raw_movie(sprintf('/Volumes/Data/2015-11-09-1/Visual/pc2015_11_09_1_fromdata009/%d.rawMovie',imov),rawMovFrames,1);
subtract_movies{3}=mean(stim,1);
subtract_movies{3}=mean(stim,1)*0+127.5;
movie=stim-repmat(subtract_movies{3},[rawMovFrames,1,1]);
movie=movie/255;

movie = permute(movie,[2,3,1]);
movie=repelem(movie,1,1,interval);


 % WN details to get the significant stixels box   
WN_datafile = '/Volumes/Analysis/2015-11-09-1/d14_36-norefit/data009/data009';
movie_xml = 'RGB-8-2-0.48-11111';
stim_length=1800;% 
% load cell IDs for ON type 

% Load spikes
datarun=load_data(WN_datafile)
datarun=load_params(datarun)
datarun = load_sta(datarun);
cellTypeId = 1;
InterestingCell_vis_id=datarun.cell_types{cellTypeId}.cell_ids; 
idata=1;
for cellID_list = InterestingCell_vis_id;
    
 % spikes for null movie
condDuration=1800;dataRuns = 25;% 25 - ON , 27 - OFF
Null_datafile = sprintf('/Volumes/Analysis/2015-11-09-1/d14_36-norefit/data0%02d',dataRuns(idata));
neuronPath = [Null_datafile,sprintf('/data0%02d.neurons',dataRuns(idata))];
[spks,spksGen_hr]=get_spikes_pc2015_11_09_1(neuronPath,cellID_list,condDuration);
   
nSU_list= [1,2,3,4,5,6,7,8,9,10];
destination= 'pc2015_11_09_1_analysis_fits/SUs_data025_exp'
sta_depth=30; % CAREFUL TO CHANGE IT TO LONG .. FOR SLOWER CELLS

[~] = get_gmlm_sta4(WN_datafile,movie_xml,stim_length,cellID_list,nSU_list,destination,sta_depth,movie,spks)
end

