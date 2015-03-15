addpath(genpath('../null_analyse/'));
addpath(genpath('../null_analyse/analyse_functions'));
startup_null_analyse_tenessee
%startup_null_analyse_bertha

%%
% Condition strings
nConditions=6;
condDuration=10.6;
cond_str=cell(6,1);
cond_str{1}='Original';
cond_str{2}='ON-OFF Mix Spatial BW';
cond_str{3}='ON-OFF Mix Spatio-Temporal BW';
cond_str{4}='ON-OFF Mix Original NSEM';
cond_str{5}='ON-OFF Mix Spatial Null NSEM ';
cond_str{6}='SBC';
interestingConditions=[1,2,3,4];

%% Load Movies
rawMovFrames=1270/(1);
figure;
icnt=0;
for imov=[1,2,4,5,6,8]
    [stim,height,width,header_size] = get_raw_movie(sprintf('/Volumes/Data/2015-03-09-7/Visual/pc2015_03_09_7_data003/%d.rawMovie',imov),rawMovFrames,1);
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
for imov=[1,2,4,5,6,8]
    [stim,height,width,header_size] = get_raw_movie(sprintf('/Volumes/Data/2015-03-09-7/Visual/pc2015_03_09_7_data003/%d.rawMovie',imov),rawMovFrames,1);
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
   hgexport(h,sprintf('/Volumes/Analysis/nishal/analyse_2015_03_09_7/data003/cMap.eps'),s);
  
%% data006 from data003


WN_datafile = '2015-03-09-7/streamed/data003/data003';
Null_datafile = '/Volumes/Analysis/2015-03-09-7/data006-from-data003_streamed_nps';
%WN_mov='/Volumes/Analysis/stimuli/white-noise-xml/RGB-10-2-0.48-11111.xml';
WN_datafile_full = '/Volumes/Analysis/2015-03-09-7/streamed/data003/data003';
neuronPath = [Null_datafile,sprintf('/data006-from-data003_streamed_nps.neurons')];
imov=14;

datarun=load_data(WN_datafile)
datarun=load_params(datarun)


% cellTypeId=[1]; % 1 for On Parasols, 2 for Off parasols
% InterestingCell_vis_id=[];
% for icellType=cellTypeId
%     icellType 
%     InterestingCell_vis_id=[InterestingCell_vis_id,datarun.cell_types{icellType}.cell_ids];
% end 

%InterestingCell_vis_id = [5568,1726,5252,3061]; % OFF 
% InterestingCell_vis_id = [6106,1835,7730]; % ON
InterestingCell_vis_id = [4501]; % SBC


NullCells1=[5568,1726,5252,3061];  % OFF
NullCells2=[6106,1835,7730];   % ON
NullCells3=[4501]; % SBC

cellType_str='SBC'
condDuration=1270/120;
nConditions=6;
for ref_cell_number=1:length(InterestingCell_vis_id); %11
    close all
    cellID=InterestingCell_vis_id(ref_cell_number);
    % [spkColl,spkCondColl,h]=plot_raster_script_pc2015_02_24_2(datarun,WN_datafile,WN_datafile_full,Null_datafile,InterestingCell_vis_id,imov,ref_cell_number,nConditions,condDuration,cond_str,neuronPath);
 [spkColl,spkCondColl,h]=plot_raster_script_pc2015_03_09_7_light(cellID,nConditions,condDuration,cond_str,neuronPath);
 
 plot_mosaic_pc2015_03_09_7(datarun,InterestingCell_vis_id,ref_cell_number,NullCells1,NullCells2,NullCells3);
   InterestingCell_vis_id(ref_cell_number)
    if(~isdir(sprintf('/Volumes/Analysis/nishal/analyse_2015_03_09_7/data006/CellType_%s',cellType_str)))
    mkdir(sprintf('/Volumes/Analysis/nishal/analyse_2015_03_09_7/data006/CellType_%s',cellType_str));
    end
   s=hgexport('readstyle','ras_mos4');
   hgexport(h,sprintf('/Volumes/Analysis/nishal/analyse_2015_03_09_7/data006/CellType_%s/CellID_%d.eps',cellType_str,InterestingCell_vis_id(ref_cell_number)),s);
  
    %testsuite_prediction
  %[timeLogData,psthData] = psth_variability(spkCondColl,nConditions,condMovies,cond_str,InterestingCell_vis_id,imov,ref_cell_number,interestingConditions);
   % pause
end

