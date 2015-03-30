addpath(genpath('../null_analyse/'));
addpath(genpath('../null_analyse/analyse_functions'));
startup_null_analyse_tenessee
%startup_null_analyse_bertha

%%
% Condition strings
nConditions=5;
condDuration=1270/120;
cond_str=cell(6,1);
cond_str{1}='Original';
cond_str{2}='Cell group 1 Spatial';
cond_str{3}='Cell group 1 Spatio-temporal';
cond_str{4}='Cell group 2 Spatial';
cond_str{5}='On Parasol';
interestingConditions=[1,2,3,4,5];

%% Load Movies
rawMovFrames=1270/(1);
figure;
icnt=0;
for imov=[1,2,4,6,8]
    [stim,height,width,header_size] = get_raw_movie(sprintf('/Volumes/Data/2015-03-09-3/Visual/Null Movies/pc2015_03_09_3_data000/%d.rawMovie',imov),rawMovFrames,1);
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

rawMovFrames=1270/(1);
figure;
icnt=0;
cMap = cell(6,1);
h=figure('Color','w');
for imov=[1,2,4,6,8]
    [stim,height,width,header_size] = get_raw_movie(sprintf('/Volumes/Data/2015-03-09-3/Visual/Null Movies/pc2015_03_09_3_data000/%d.rawMovie',imov),rawMovFrames,1);
    subtract_movies{3}=mean(stim,1);
    subtract_movies{3}=mean(stim,1)*0+127.5;
    movie=stim-repmat(subtract_movies{3},[rawMovFrames,1,1]);
    movie=movie/255;

    icnt=icnt+1;
    
    subplot(3,2,icnt);
    qq=movie;
   
    cMap{icnt}=contrastMap(qq);
   
    imagesc(cMap{icnt});
    %caxis([3,6]);
    colorbar
    axis image
    title(sprintf('cMap: %d',imov));
end

   s=hgexport('readstyle','cMap');
   hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_3/data000/cMap.eps'),s);
  
%% data003 from data000


WN_datafile = '2015-03-09-3/streamed/data000/data000';
Null_datafile = '/Volumes/Analysis/2015-03-09-3/data003-04-from-data000_streamed_nps';
%WN_mov='/Volumes/Analysis/stimuli/white-noise-xml/RGB-10-2-0.48-11111.xml';
WN_datafile_full = '/Volumes/Analysis/2015-03-09-3/streamed/data000/data000';
neuronPath = [Null_datafile,sprintf('/data003-04-from-data000_streamed_nps.neurons')];


datarun=load_data(WN_datafile)
datarun=load_params(datarun)


cellTypeId=[1]; % 1 for On Parasols, 2 for Off parasols
InterestingCell_vis_id=[];
for icellType=cellTypeId
    icellType 
    InterestingCell_vis_id=[InterestingCell_vis_id,datarun.cell_types{icellType}.cell_ids];
end 

%InterestingCell_vis_id = [5568,1726,5252,3061]; % OFF 
% InterestingCell_vis_id = [6106,1835,7730]; % ON
%InterestingCell_vis_id = [4501]; % SBC


NullCells1=[7368,872,4114,2448,5581];  % OFF
NullCells2=[4863,2615,3783,872,5735];   % ON


condDuration=1270/120;
nConditions=5;
for ref_cell_number=1:length(InterestingCell_vis_id); %11
    close all
    cellID=InterestingCell_vis_id(ref_cell_number);
    % [spkColl,spkCondColl,h]=plot_raster_script_pc2015_02_24_2(datarun,WN_datafile,WN_datafile_full,Null_datafile,InterestingCell_vis_id,imov,ref_cell_number,nConditions,condDuration,cond_str,neuronPath);
 [spkColl,spkCondColl,h]=plot_raster_script_pc2015_03_09_3_light(cellID,nConditions,condDuration,cond_str,neuronPath);
 
 plot_mosaic_pc2015_03_09_3(datarun,InterestingCell_vis_id,ref_cell_number,NullCells1,NullCells2);
   InterestingCell_vis_id(ref_cell_number)
    if(~isdir(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_3/data003/CellType_%s',datarun.cell_types{cellTypeId}.name)))
    mkdir(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_3/data003/CellType_%s',datarun.cell_types{cellTypeId}.name));
    end
   s=hgexport('readstyle','ras_mos4');
   hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_3/data003/CellType_%s/CellID_%d.eps',datarun.cell_types{cellTypeId}.name,InterestingCell_vis_id(ref_cell_number)),s);
  
    %testsuite_prediction
  %[timeLogData,psthData] = psth_variability(spkCondColl,nConditions,condMovies,cond_str,InterestingCell_vis_id,imov,ref_cell_number,interestingConditions);
   % pause
end



%% analyze data012,data013,data008
% Condition strings
nConditions=3;
condDuration=1272/120;
cond_str=cell(3,1);
cond_str{1}='Original';
cond_str{2}='Cell group 1 Spatial';
cond_str{3}='OFF parasol';
interestingConditions=[1,2,3];

%% data012-from-data008
WN_datafile = '2015-03-09-3/streamed/data008/data008';
Null_datafile = '/Volumes/Analysis/2015-03-09-3/data012-13-from-data008_streamed_nps';
%WN_mov='/Volumes/Analysis/stimuli/white-noise-xml/RGB-10-2-0.48-11111.xml';
WN_datafile_full = '/Volumes/Analysis/2015-03-09-3/streamed/data008/data008';
neuronPath = [Null_datafile,sprintf('/data012-13-from-data008_streamed_nps.neurons')];


datarun=load_data(WN_datafile)
datarun=load_params(datarun)


cellTypeId=[2]; % 1 for On Parasols, 2 for Off parasols
InterestingCell_vis_id=[];
for icellType=cellTypeId
    icellType 
    InterestingCell_vis_id=[InterestingCell_vis_id,datarun.cell_types{icellType}.cell_ids];
end 

InterestingCell_vis_id = [4083,1804,2448,5221,1068]; % OFF 


NullCells1=[4083,1804,2448,5221,1068];  % OFF
NullCells2=[];   % ON


condDuration=1272/120;
nConditions=3;
for ref_cell_number=1:length(InterestingCell_vis_id); %11
    close all
    cellID=InterestingCell_vis_id(ref_cell_number);
    % [spkColl,spkCondColl,h]=plot_raster_script_pc2015_02_24_2(datarun,WN_datafile,WN_datafile_full,Null_datafile,InterestingCell_vis_id,imov,ref_cell_number,nConditions,condDuration,cond_str,neuronPath);
 [spkColl,spkCondColl,h]=plot_raster_script_pc2015_03_09_3_light(cellID,nConditions,condDuration,cond_str,neuronPath);
 
 plot_mosaic_pc2015_03_09_3(datarun,InterestingCell_vis_id,ref_cell_number,NullCells1,NullCells2);
   InterestingCell_vis_id(ref_cell_number)
    if(~isdir(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_3/data012/CellType_%s/CellID_%d',datarun.cell_types{cellTypeId}.name,cellID)))
    mkdir(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_3/data012/CellType_%s/CellID_%d',datarun.cell_types{cellTypeId}.name,cellID));
    end
   s=hgexport('readstyle','ras_mos4');
   hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_3/data012/CellType_%s/CellID_%d/recorded.eps',datarun.cell_types{cellTypeId}.name,InterestingCell_vis_id(ref_cell_number)),s);
  
    %testsuite_prediction
  %[timeLogData,psthData] = psth_variability(spkCondColl,nConditions,condMovies,cond_str,InterestingCell_vis_id,imov,ref_cell_number,interestingConditions);
   % pause
end

