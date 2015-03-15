addpath(genpath('../null_analyse/'));
addpath(genpath('../null_analyse/analyse_functions'));
startup_null_analyse_tenessee
%startup_null_analyse_bertha

%%
% Condition strings
nConditions=3;
condDuration=12;
cond_str=cell(6,1);
cond_str{1}='Original';
cond_str{2}='On selected ';
cond_str{3}='All ON';
cond_str{4}='Off selected ';
cond_str{5}='All OFF';
cond_str{6}='SBC';

interestingConditions=[1,2,3,4,5,6];
%% Load Movies
rawMovFrames=1272/(6);
figure;
icnt=0;
for imov=[1,2,4,6,8,10]
    [stim,height,width,header_size] = get_raw_movie(sprintf('/Volumes/Data/2015-03-09-2/Visual/pc2015_03_09_2_data038/%d.rawMovie',imov),rawMovFrames,1);
    subtract_movies{3}=mean(stim,1);
    subtract_movies{3}=mean(stim,1)*0+127.5;
    movie=stim-repmat(subtract_movies{3},[rawMovFrames,1,1]);
    movie=movie/255;

    icnt=icnt+1;
    
    subplot(3,2,icnt);
    qq=movie;
    hist(qq(:),20)
    title(sprintf('Movie pixel histogram %d',imov));
end

rawMovFrames=1272/(6);
figure;
icnt=0;
cMap = cell(6,1);
h=figure('Color','w');
for imov=[1,2,4,6,8,10]
    [stim,height,width,header_size] = get_raw_movie(sprintf('/Volumes/Data/2015-03-09-2/Visual/pc2015_03_09_2_data038/%d.rawMovie',imov),rawMovFrames,1);
    subtract_movies{3}=mean(stim,1);
    subtract_movies{3}=mean(stim,1)*0+127.5;
    movie=stim-repmat(subtract_movies{3},[rawMovFrames,1,1]);
    movie=movie/255;

    icnt=icnt+1;
    
    subplot(3,2,icnt);
    qq=movie;
   
    cMap{icnt}=contrastMap(qq);
   
    imagesc(cMap{icnt});
    caxis([3,6]);
    colorbar
    axis image
    title(sprintf('cMap: %d',imov));
end

   s=hgexport('readstyle','cMap');
   hgexport(h,sprintf('/Volumes/Analysis/nishal/analyse_2015_03_09_2/data038/cMap.eps'),s);
  
%% data041 from data031


WN_datafile = '2015-03-09-2/data031/data031';
Null_datafile = '/Volumes/Analysis/2015-03-09-2/data041-from-data031_nps';
%WN_mov='/Volumes/Analysis/stimuli/white-noise-xml/RGB-10-2-0.48-11111.xml';
WN_datafile_full = '/Volumes/Analysis/2015-03-09-2/data031/data031';
neuronPath = [Null_datafile,sprintf('/data041-from-data031_nps.neurons')];
imov=14;

datarun=load_data(WN_datafile)
datarun=load_params(datarun)


cellTypeId=[1]; % 1 for On Parasols, 2 for Off parasols
InterestingCell_vis_id=[];
for icellType=cellTypeId
    icellType 
    InterestingCell_vis_id=[InterestingCell_vis_id,datarun.cell_types{icellType}.cell_ids];
end 

cellTypeUsed=cellTypeId*ones(length(InterestingCell_vis_id),1);
NullCells1=[2206,4773,1008]; % cell : 6826 ? 
NullCells2=datarun.cell_types{1}.cell_ids;
NullCells3=[2461,1021,7188,767,4502,5911];
NullCells4=datarun.cell_types{2}.cell_ids;
NullCells5=[5732];


condDuration=10.6;
nConditions=6;
for ref_cell_number=1:length(InterestingCell_vis_id); %11
    close all
    cellID=InterestingCell_vis_id(ref_cell_number);
    % [spkColl,spkCondColl,h]=plot_raster_script_pc2015_02_24_2(datarun,WN_datafile,WN_datafile_full,Null_datafile,InterestingCell_vis_id,imov,ref_cell_number,nConditions,condDuration,cond_str,neuronPath);
 [spkColl,spkCondColl,h]=plot_raster_script_pc2015_03_09_2_light(cellID,nConditions,condDuration,cond_str,neuronPath);
 
 plot_mosaic_pc2015_03_09_2(datarun,InterestingCell_vis_id,ref_cell_number,NullCells1,NullCells2);
   InterestingCell_vis_id(ref_cell_number)
    if(~isdir(sprintf('/Volumes/Analysis/nishal/analyse_2015_03_09_2/data041/CellType_%s',datarun.cell_types{cellTypeId}.name)))
    mkdir(sprintf('/Volumes/Analysis/nishal/analyse_2015_03_09_2/data041/CellType_%s',datarun.cell_types{cellTypeId}.name));
    end
   s=hgexport('readstyle','ras_mos4');
   hgexport(h,sprintf('/Volumes/Analysis/nishal/analyse_2015_03_09_2/data041/CellType_%s/CellID_%d.eps',datarun.cell_types{cellTypeId}.name,InterestingCell_vis_id(ref_cell_number)),s);
  
    %testsuite_prediction
  %[timeLogData,psthData] = psth_variability(spkCondColl,nConditions,condMovies,cond_str,InterestingCell_vis_id,imov,ref_cell_number,interestingConditions);
   % pause
end


%% data042 from data038


WN_datafile = '2015-03-09-2/streamed/data038/data038';
Null_datafile = '/Volumes/Analysis/2015-03-09-2/data042-from-data038_streamed_nps';
%WN_mov='/Volumes/Analysis/stimuli/white-noise-xml/RGB-10-2-0.48-11111.xml';
WN_datafile_full = '/Volumes/Analysis/2015-03-09-2/data038/data038';
neuronPath = [Null_datafile,sprintf('/data042-from-data038_streamed_nps.neurons')];


datarun=load_data(WN_datafile)
datarun=load_params(datarun)


cellTypeId=[1]; % 1 for On Parasols, 2 for Off parasols
InterestingCell_vis_id=[];
for icellType=cellTypeId
    icellType 
    InterestingCell_vis_id=[InterestingCell_vis_id,datarun.cell_types{icellType}.cell_ids];
end 

cellTypeUsed=cellTypeId*ones(length(InterestingCell_vis_id),1);
NullCells1=[2206,4548,6826,662]; % cell : 6826 ? 
NullCells2=datarun.cell_types{1}.cell_ids;
NullCells3=[2596,1232,7172,842,4547,5911];
NullCells4=datarun.cell_types{2}.cell_ids;
NullCells5=[5768];


condDuration=10.6;
nConditions=6;
for ref_cell_number=1:length(InterestingCell_vis_id); %11
    close all
    cellID=InterestingCell_vis_id(ref_cell_number);
    % [spkColl,spkCondColl,h]=plot_raster_script_pc2015_02_24_2(datarun,WN_datafile,WN_datafile_full,Null_datafile,InterestingCell_vis_id,imov,ref_cell_number,nConditions,condDuration,cond_str,neuronPath);
 [spkColl,spkCondColl,h]=plot_raster_script_pc2015_03_09_2_light(cellID,nConditions,condDuration,cond_str,neuronPath);
 
 plot_mosaic_pc2015_03_09_2(datarun,InterestingCell_vis_id,ref_cell_number,NullCells1,NullCells2);
   InterestingCell_vis_id(ref_cell_number)
    if(~isdir(sprintf('/Volumes/Analysis/nishal/analyse_2015_03_09_2/data042/CellType_%s',datarun.cell_types{cellTypeId}.name)))
    mkdir(sprintf('/Volumes/Analysis/nishal/analyse_2015_03_09_2/data042/CellType_%s',datarun.cell_types{cellTypeId}.name));
    end
   s=hgexport('readstyle','ras_mos4');
   hgexport(h,sprintf('/Volumes/Analysis/nishal/analyse_2015_03_09_2/data042/CellType_%s/CellID_%d.eps',datarun.cell_types{cellTypeId}.name,InterestingCell_vis_id(ref_cell_number)),s);
  
    %testsuite_prediction
  %[timeLogData,psthData] = psth_variability(spkCondColl,nConditions,condMovies,cond_str,InterestingCell_vis_id,imov,ref_cell_number,interestingConditions);
   % pause
end

%