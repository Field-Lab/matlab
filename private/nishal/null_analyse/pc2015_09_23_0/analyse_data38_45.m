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

rawMovFrames=1200/(4);
figure;
icnt=0;
% make pixel histogram
for imov=[1,2,4,10,12,5,6,8]
    [stim,height,width,header_size] = get_raw_movie(sprintf('/Volumes/Data/2015-09-23-0/Visual/2015-09-23-0_data017/%d.rawMovie',imov),rawMovFrames,1);
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

% make movies
interval=4;
condMov=cell(nConditions,1);
rawMovFrames=1200/(4);
icnt=0;
% make pixel histogram
for imov=[1,2,4,10,12,5,6,8]
    [stim,height,width,header_size] = get_raw_movie(sprintf('/Volumes/Data/2015-09-23-0/Visual/2015-09-23-0_data017/%d.rawMovie',imov),rawMovFrames,1);
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
    
end

% make contrast map
rawMovFrames=1200/(4);
figure;
icnt=0;
cMap = cell(8,1);
h=figure('Color','w');
for imov=[1,2,4,10,12,5,6,8]
    [stim,height,width,header_size] = get_raw_movie(sprintf('/Volumes/Data/2015-09-23-0/Visual/2015-09-23-0_data017/%d.rawMovie',imov),rawMovFrames,1);
    subtract_movies{3}=mean(stim,1);
    subtract_movies{3}=mean(stim,1)*0+127.5;
    movie=stim-repmat(subtract_movies{3},[rawMovFrames,1,1]);
    movie=movie/255;
    
    icnt=icnt+1;
    
    subplot(34,2,icnt);
    qq=movie;
    
    cMap{icnt}=contrastMap(qq);
    
    imagesc(cMap{icnt});
    caxis([3,6]);
    colorbar
    axis image
    title(sprintf('cMap: %d',imov));
end

s=hgexport('readstyle','cMap');
hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_09_23_0/data000/cMap.eps'),s);


%% data041 from data031


WN_datafile = '/Volumes/Analysis/2015-09-23-0/d17_25_46-norefit/data017-from-data017_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035_data036_data037_data038_data039_data040_data041_data042_data043_data044_data045_data046/data017-from-data017_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035_data036_data037_data038_data039_data040_data041_data042_data043_data044_data045_data046';
%WN_mov='/Volumes/Analysis/stimuli/white-noise-xml/RGB-10-2-0.48-11111.xml';
imov=14;
WN_datafile_full = '/Volumes/Analysis/2015-09-23-0/streamed/data017/';

datarun=load_data(WN_datafile)
datarun=load_params(datarun)

cellTypeId = 2;
InterestingCell_vis_id=datarun.cell_types{cellTypeId}.cell_ids; %[556,1278,1384,407,1516,2150,2401,3361,4066,4683,5611,6106,6005,7246,7562,3946];
cellTypeUsed=cellTypeId*ones(length(InterestingCell_vis_id),1);

condDuration=10;
nConditions=1;

cols='rkrkrkrkrkrkkrkrkrkr';
spkCondColl=cell(8,1);

for ref_cell_number=1:length(InterestingCell_vis_id); %11
    close all
    
    %data038
    Null_datafile = '/Volumes/Analysis/2015-09-23-0/d17_25_46-norefit/data038-from-data017_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035_data036_data037_data038_data039_data040_data041_data042_data043_data044_data045_data046';
    neuronPath = [Null_datafile,sprintf('/data038-from-data017_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035_data036_data037_data038_data039_data040_data041_data042_data043_data044_data045_data046.neurons')];
    cellID=InterestingCell_vis_id(ref_cell_number);
    [spkColl,spkCondColl{1},h]=plot_raster_script_pc2015_09_23_0_light(cellID,nConditions,condDuration,cond_str,neuronPath);
   
      %data039
    Null_datafile = '/Volumes/Analysis/2015-09-23-0/d17_25_46-norefit/data039-from-data017_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035_data036_data037_data038_data039_data040_data041_data042_data043_data044_data045_data046';
    neuronPath = [Null_datafile,sprintf('/data039-from-data017_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035_data036_data037_data038_data039_data040_data041_data042_data043_data044_data045_data046.neurons')];
  cellID=InterestingCell_vis_id(ref_cell_number);
    [spkColl,spkCondColl{2},h]=plot_raster_script_pc2015_09_23_0_light(cellID,nConditions,condDuration,cond_str,neuronPath);
   
    
       %data040
    Null_datafile = '/Volumes/Analysis/2015-09-23-0/d17_25_46-norefit/data040-from-data017_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035_data036_data037_data038_data039_data040_data041_data042_data043_data044_data045_data046';
    neuronPath = [Null_datafile,sprintf('/data040-from-data017_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035_data036_data037_data038_data039_data040_data041_data042_data043_data044_data045_data046.neurons')];
   cellID=InterestingCell_vis_id(ref_cell_number);
    [spkColl,spkCondColl{3},h]=plot_raster_script_pc2015_09_23_0_light(cellID,nConditions,condDuration,cond_str,neuronPath);
   
    
    
        %data041
    Null_datafile = '/Volumes/Analysis/2015-09-23-0/d17_25_46-norefit/data041-from-data017_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035_data036_data037_data038_data039_data040_data041_data042_data043_data044_data045_data046';
    neuronPath = [Null_datafile,sprintf('/data041-from-data017_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035_data036_data037_data038_data039_data040_data041_data042_data043_data044_data045_data046.neurons')];
  cellID=InterestingCell_vis_id(ref_cell_number);
    [spkColl,spkCondColl{4},h]=plot_raster_script_pc2015_09_23_0_light(cellID,nConditions,condDuration,cond_str,neuronPath);
   
    
    
     %data042
 Null_datafile = '/Volumes/Analysis/2015-09-23-0/d17_25_46-norefit/data042-from-data017_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035_data036_data037_data038_data039_data040_data041_data042_data043_data044_data045_data046';
    neuronPath = [Null_datafile,sprintf('/data042-from-data017_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035_data036_data037_data038_data039_data040_data041_data042_data043_data044_data045_data046.neurons')];
     cellID=InterestingCell_vis_id(ref_cell_number);
    [spkColl,spkCondColl{5},h]=plot_raster_script_pc2015_09_23_0_light(cellID,nConditions,condDuration,cond_str,neuronPath);
   
     %data043
    Null_datafile = '/Volumes/Analysis/2015-09-23-0/d17_25_46-norefit/data043-from-data017_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035_data036_data037_data038_data039_data040_data041_data042_data043_data044_data045_data046';
    neuronPath = [Null_datafile,sprintf('/data043-from-data017_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035_data036_data037_data038_data039_data040_data041_data042_data043_data044_data045_data046.neurons')];
     cellID=InterestingCell_vis_id(ref_cell_number);
    [spkColl,spkCondColl{6},h]=plot_raster_script_pc2015_09_23_0_light(cellID,nConditions,condDuration,cond_str,neuronPath);
   
    
    %data044
    Null_datafile = '/Volumes/Analysis/2015-09-23-0/d17_25_46-norefit/data044-from-data017_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035_data036_data037_data038_data039_data040_data041_data042_data043_data044_data045_data046';
    neuronPath = [Null_datafile,sprintf('/data044-from-data017_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035_data036_data037_data038_data039_data040_data041_data042_data043_data044_data045_data046.neurons')];
      cellID=InterestingCell_vis_id(ref_cell_number);
    [spkColl,spkCondColl{7},h]=plot_raster_script_pc2015_09_23_0_light(cellID,nConditions,condDuration,cond_str,neuronPath);
   
      %data045
    Null_datafile = '/Volumes/Analysis/2015-09-23-0/d17_25_46-norefit/data045-from-data017_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035_data036_data037_data038_data039_data040_data041_data042_data043_data044_data045_data046';
    neuronPath = [Null_datafile,sprintf('/data045-from-data017_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035_data036_data037_data038_data039_data040_data041_data042_data043_data044_data045_data046.neurons')];
    cellID=InterestingCell_vis_id(ref_cell_number);
    [spkColl,spkCondColl{8},h]=plot_raster_script_pc2015_09_23_0_light(cellID,nConditions,condDuration,cond_str,neuronPath);
 
    
    
    h=figure;
    for irun = 1:8
    plot(spkCondColl{irun}.xPoints/20000,spkCondColl{irun}.yPoints - (irun-1)*30,cols(irun));
    hold on;
    end
    set(gca,'yTick',[]);
    ylim([-7*30,30]);
    InterestingCell_vis_id(ref_cell_number)
    
    if(~isdir(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_09_23_0/d38_45/CellType_%s',datarun.cell_types{cellTypeId}.name)))
        mkdir(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_09_23_0/d38_45/CellType_%s',datarun.cell_types{cellTypeId}.name));
    end
   s=hgexport('readstyle','raster');
   hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_09_23_0/d38_45/CellType_%s/CellID_%d.eps',datarun.cell_types{cellTypeId}.name,InterestingCell_vis_id(ref_cell_number)),s);
end


