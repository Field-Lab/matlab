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

interestingConditions=[1,2,3,4,5,6,7];
%% Load Movies
% 
rawMovFrames=1200/(4);
figure;
icnt=0;
% make pixel histogram
for imov=[16,18,19,21,23,24,26]
    [stim,height,width,header_size] = get_raw_movie(sprintf('/Volumes/Data/2015-10-06-0/Visual/null/%d.rawMovie',imov),rawMovFrames,1);
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
% % make contrast map
% rawMovFrames=1200/(4);
% figure;
% icnt=0;
% cMap = cell(8,1);
% h=figure('Color','w');
% for imov=[1,2,4,10,12,5,6,8]
%     [stim,height,width,header_size] = get_raw_movie(sprintf('/Volumes/Data/2015-09-23-0/Visual/2015-09-23-0_data017/%d.rawMovie',imov),rawMovFrames,1);
%     subtract_movies{3}=mean(stim,1);
%     subtract_movies{3}=mean(stim,1)*0+127.5;
%     movie=stim-repmat(subtract_movies{3},[rawMovFrames,1,1]);
%     movie=movie/255;
%     
%     icnt=icnt+1;
%     
%     subplot(34,2,icnt);
%     qq=movie;
%     
%     cMap{icnt}=contrastMap(qq);
%     
%     imagesc(cMap{icnt});
%     caxis([3,6]);
%     colorbar
%     axis image
%     title(sprintf('cMap: %d',imov));
% end
% 
% s=hgexport('readstyle','cMap');
% hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_09_23_0/data000/cMap.eps'),s);


%% correct SU nulled - data016,18,17,28,30,32,34,27

WN_datafile = '/Volumes/Analysis/2015-10-06-0/d14_35-norefit/data014-from-data014_data016_data017_data018_data019_data020_data021_data022_data023_data024_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035/data014-from-data014_data016_data017_data018_data019_data020_data021_data022_data023_data024_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035';
%WN_mov='/Volumes/Analysis/stimuli/white-noise-xml/RGB-10-2-0.48-11111.xml';
imov=14;
WN_datafile_full = '/Volumes/Analysis/2015-10-06-0/streamed/data014/';

datarun=load_data(WN_datafile)
datarun=load_params(datarun)

cellTypeId = 8%8,1;
InterestingCell_vis_id=datarun.cell_types{cellTypeId}.cell_ids; %[556,1278,1384,407,1516,2150,2401,3361,4066,4683,5611,6106,6005,7246,7562,3946];
cellTypeUsed=cellTypeId*ones(length(InterestingCell_vis_id),1);

condDuration=10;
nConditions=1;

cols='rkrkrkrkrkrkkrkrkrkr';
spkCondColl=cell(8,1);

for ref_cell_number=1:length(InterestingCell_vis_id); %11
    close all
    cellID=InterestingCell_vis_id(ref_cell_number)
        
    %data016
    Null_datafile = '/Volumes/Analysis/2015-10-06-0/d14_35-norefit/data016-from-data014_data016_data017_data018_data019_data020_data021_data022_data023_data024_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035';
    neuronPath = [Null_datafile,sprintf('/data016-from-data014_data016_data017_data018_data019_data020_data021_data022_data023_data024_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035.neurons')];
    [spkColl,spkCondColl{1},h]=plot_raster_script_pc2015_09_23_0_light(cellID,nConditions,condDuration,cond_str,neuronPath);
   
    
    %data018
    Null_datafile = '/Volumes/Analysis/2015-10-06-0/d14_35-norefit/data018-from-data014_data016_data017_data018_data019_data020_data021_data022_data023_data024_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035';
    neuronPath = [Null_datafile,sprintf('/data018-from-data014_data016_data017_data018_data019_data020_data021_data022_data023_data024_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035.neurons')];
    [spkColl,spkCondColl{2},h]=plot_raster_script_pc2015_09_23_0_light(cellID,nConditions,condDuration,cond_str,neuronPath);
    
      %data017
    Null_datafile = '/Volumes/Analysis/2015-10-06-0/d14_35-norefit/data017-from-data014_data016_data017_data018_data019_data020_data021_data022_data023_data024_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035';
    neuronPath = [Null_datafile,sprintf('/data017-from-data014_data016_data017_data018_data019_data020_data021_data022_data023_data024_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035.neurons')];
    [spkColl,spkCondColl{3},h]=plot_raster_script_pc2015_09_23_0_light(cellID,nConditions,condDuration,cond_str,neuronPath);
   
    
    
    %data028
    Null_datafile = '/Volumes/Analysis/2015-10-06-0/d14_35-norefit/data028-from-data014_data016_data017_data018_data019_data020_data021_data022_data023_data024_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035';
    neuronPath = [Null_datafile,sprintf('/data028-from-data014_data016_data017_data018_data019_data020_data021_data022_data023_data024_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035.neurons')];
    [spkColl,spkCondColl{4},h]=plot_raster_script_pc2015_09_23_0_light(cellID,nConditions,condDuration,cond_str,neuronPath);
   
    %data030
    Null_datafile = '/Volumes/Analysis/2015-10-06-0/d14_35-norefit/data030-from-data014_data016_data017_data018_data019_data020_data021_data022_data023_data024_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035';
    neuronPath = [Null_datafile,sprintf('/data030-from-data014_data016_data017_data018_data019_data020_data021_data022_data023_data024_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035.neurons')];
   [spkColl,spkCondColl{5},h]=plot_raster_script_pc2015_09_23_0_light(cellID,nConditions,condDuration,cond_str,neuronPath);
   
    
    %data032
    Null_datafile = '/Volumes/Analysis/2015-10-06-0/d14_35-norefit/data032-from-data014_data016_data017_data018_data019_data020_data021_data022_data023_data024_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035';
    neuronPath = [Null_datafile,sprintf('/data032-from-data014_data016_data017_data018_data019_data020_data021_data022_data023_data024_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035.neurons')];
   [spkColl,spkCondColl{6},h]=plot_raster_script_pc2015_09_23_0_light(cellID,nConditions,condDuration,cond_str,neuronPath);
   
    %data034
    Null_datafile = '/Volumes/Analysis/2015-10-06-0/d14_35-norefit/data034-from-data014_data016_data017_data018_data019_data020_data021_data022_data023_data024_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035';
    neuronPath = [Null_datafile,sprintf('/data034-from-data014_data016_data017_data018_data019_data020_data021_data022_data023_data024_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035.neurons')];
    [spkColl,spkCondColl{7},h]=plot_raster_script_pc2015_09_23_0_light(cellID,nConditions,condDuration,cond_str,neuronPath);
 
    
   
    
    
     %data027
    Null_datafile = '/Volumes/Analysis/2015-10-06-0/d14_35-norefit/data027-from-data014_data016_data017_data018_data019_data020_data021_data022_data023_data024_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035';
    neuronPath = [Null_datafile,sprintf('/data027-from-data014_data016_data017_data018_data019_data020_data021_data022_data023_data024_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035.neurons')];
    [spkColl,spkCondColl{8},h]=plot_raster_script_pc2015_09_23_0_light(cellID,nConditions,condDuration,cond_str,neuronPath);
   
    
    
    h=figure;
    for irun = 1:8
    plot(spkCondColl{irun}.xPoints/20000,spkCondColl{irun}.yPoints - (irun-1)*30,cols(irun));
    hold on;
    end
    set(gca,'yTick',[]);
    ylim([-7*30,30]);
    InterestingCell_vis_id(ref_cell_number)
    
    if(~isdir(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_10_06_0/dsu_null_28_30_32_34/CellType_%s',datarun.cell_types{cellTypeId}.name)))
        mkdir(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_10_06_0/dsu_null_28_30_32_34/CellType_%s',datarun.cell_types{cellTypeId}.name));
    end
   s=hgexport('readstyle','raster');
   hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_10_06_0/dsu_null_28_30_32_34/CellType_%s/CellID_%d.eps',datarun.cell_types{cellTypeId}.name,InterestingCell_vis_id(ref_cell_number)),s);
end

%% wrong SU nulled - data016,18,17,29,31,33,35,27

dataRuns = [16,17,18,29,31,33,35,27];

WN_datafile = '/Volumes/Analysis/2015-10-06-0/d14_35-norefit/data014-from-data014_data016_data017_data018_data019_data020_data021_data022_data023_data024_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035/data014-from-data014_data016_data017_data018_data019_data020_data021_data022_data023_data024_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035';
%WN_mov='/Volumes/Analysis/stimuli/white-noise-xml/RGB-10-2-0.48-11111.xml';
imov=14;
WN_datafile_full = '/Volumes/Analysis/2015-10-06-0/streamed/data014/';

datarun=load_data(WN_datafile)
datarun=load_params(datarun)

cellTypeId = 1%8,1;
InterestingCell_vis_id=datarun.cell_types{cellTypeId}.cell_ids; %[556,1278,1384,407,1516,2150,2401,3361,4066,4683,5611,6106,6005,7246,7562,3946];
cellTypeUsed=cellTypeId*ones(length(InterestingCell_vis_id),1);

condDuration=10;
nConditions=1;

cols='rkrkrkrkrkrkkrkrkrkr';
spkCondColl=cell(8,1);

for ref_cell_number=1:length(InterestingCell_vis_id); %11
    close all
    cellID=InterestingCell_vis_id(ref_cell_number)
        
    for idata=1:length(dataRuns)
    Null_datafile = sprintf('/Volumes/Analysis/2015-10-06-0/d14_35-norefit/data0%d-from-data014_data016_data017_data018_data019_data020_data021_data022_data023_data024_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035',dataRuns(idata));
    neuronPath = [Null_datafile,sprintf('/data0%d-from-data014_data016_data017_data018_data019_data020_data021_data022_data023_data024_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035.neurons',dataRuns(idata))];
    [spkColl,spkCondColl{idata},h]=plot_raster_script_pc2015_09_23_0_light(cellID,nConditions,condDuration,cond_str,neuronPath);
    end
    
    
    h=figure;
    for irun = 1:8
    plot(spkCondColl{irun}.xPoints/20000,spkCondColl{irun}.yPoints - (irun-1)*30,cols(irun));
    hold on;
    end
    set(gca,'yTick',[]);
    ylim([-7*30,30]);
    InterestingCell_vis_id(ref_cell_number)
    
    if(~isdir(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_10_06_0/dsu_null_29_31_33_35/CellType_%s',datarun.cell_types{cellTypeId}.name)))
        mkdir(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_10_06_0/dsu_null_29_31_33_35/CellType_%s',datarun.cell_types{cellTypeId}.name));
    end
   s=hgexport('readstyle','raster');
   hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_10_06_0/dsu_null_29_31_33_35/CellType_%s/CellID_%d.eps',datarun.cell_types{cellTypeId}.name,InterestingCell_vis_id(ref_cell_number)),s);
end


%% additivity experiment - data016,17,19,20,22,24,25 - OFF

dataRuns = [16,17,19,20,22,24,25];

WN_datafile = '/Volumes/Analysis/2015-10-06-0/d14_35-norefit/data014-from-data014_data016_data017_data018_data019_data020_data021_data022_data023_data024_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035/data014-from-data014_data016_data017_data018_data019_data020_data021_data022_data023_data024_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035';
%WN_mov='/Volumes/Analysis/stimuli/white-noise-xml/RGB-10-2-0.48-11111.xml';
imov=14;
WN_datafile_full = '/Volumes/Analysis/2015-10-06-0/streamed/data014/';

datarun=load_data(WN_datafile)
datarun=load_params(datarun)

cellTypeId = 1%8,1;
InterestingCell_vis_id=datarun.cell_types{cellTypeId}.cell_ids; %[556,1278,1384,407,1516,2150,2401,3361,4066,4683,5611,6106,6005,7246,7562,3946];
cellTypeUsed=cellTypeId*ones(length(InterestingCell_vis_id),1);

condDuration=10;
nConditions=1;

cols='rkrkrkrkrkrkkrkrkrkr';
spkCondColl=cell(length(dataRuns),1);

for ref_cell_number=1:length(InterestingCell_vis_id); %11
    close all
    cellID=InterestingCell_vis_id(ref_cell_number)
        
    for idata=1:length(dataRuns)
    Null_datafile = sprintf('/Volumes/Analysis/2015-10-06-0/d14_35-norefit/data0%d-from-data014_data016_data017_data018_data019_data020_data021_data022_data023_data024_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035',dataRuns(idata));
    neuronPath = [Null_datafile,sprintf('/data0%d-from-data014_data016_data017_data018_data019_data020_data021_data022_data023_data024_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035.neurons',dataRuns(idata))];
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
    
    if(~isdir(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_10_06_0/dsu_null_16,17,19,20,22,24,25/CellType_%s',datarun.cell_types{cellTypeId}.name)))
        mkdir(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_10_06_0/dsu_null_16,17,19,20,22,24,25/CellType_%s',datarun.cell_types{cellTypeId}.name));
    end
   s=hgexport('readstyle','raster');
   hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_10_06_0/dsu_null_16,17,19,20,22,24,25/CellType_%s/CellID_%d.eps',datarun.cell_types{cellTypeId}.name,InterestingCell_vis_id(ref_cell_number)),s);
end

%% additivity experiment - data016,18,19,21,23,24,26 - ON

dataRuns = [16,18,19,21,23,24,26];

WN_datafile = '/Volumes/Analysis/2015-10-06-0/d14_35-norefit/data014-from-data014_data016_data017_data018_data019_data020_data021_data022_data023_data024_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035/data014-from-data014_data016_data017_data018_data019_data020_data021_data022_data023_data024_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035';
%WN_mov='/Volumes/Analysis/stimuli/white-noise-xml/RGB-10-2-0.48-11111.xml';
imov=14;
WN_datafile_full = '/Volumes/Analysis/2015-10-06-0/streamed/data014/';

datarun=load_data(WN_datafile)
datarun=load_params(datarun)

cellTypeId = 1%8,1;
InterestingCell_vis_id=datarun.cell_types{cellTypeId}.cell_ids; %[556,1278,1384,407,1516,2150,2401,3361,4066,4683,5611,6106,6005,7246,7562,3946];
cellTypeUsed=cellTypeId*ones(length(InterestingCell_vis_id),1);

condDuration=10;
nConditions=1;

cols='rkrkrkrkrkrkkrkrkrkr';
spkCondColl=cell(length(dataRuns),1);

for ref_cell_number=1:length(InterestingCell_vis_id); %11
    close all
    cellID=InterestingCell_vis_id(ref_cell_number)
        
    for idata=1:length(dataRuns)
    Null_datafile = sprintf('/Volumes/Analysis/2015-10-06-0/d14_35-norefit/data0%d-from-data014_data016_data017_data018_data019_data020_data021_data022_data023_data024_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035',dataRuns(idata));
    neuronPath = [Null_datafile,sprintf('/data0%d-from-data014_data016_data017_data018_data019_data020_data021_data022_data023_data024_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035.neurons',dataRuns(idata))];
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
    
    if(~isdir(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_10_06_0/dsu_null_16,18,19,21,23,24,26/CellType_%s',datarun.cell_types{cellTypeId}.name)))
        mkdir(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_10_06_0/dsu_null_16,18,19,21,23,24,26/CellType_%s',datarun.cell_types{cellTypeId}.name));
    end
   s=hgexport('readstyle','raster');
   hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_10_06_0/dsu_null_16,18,19,21,23,24,26/CellType_%s/CellID_%d.eps',datarun.cell_types{cellTypeId}.name,InterestingCell_vis_id(ref_cell_number)),s);
end










