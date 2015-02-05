addpath(genpath('../null_analyse/'));
addpath(genpath('../null_analyse/analyse_functions'));
%startup_null_analyse_tenessee
startup_null_analyse_bertha

%%
% Condition strings
nConditions=2;
condDuration=4;
cond_str=cell(2,1);
cond_str{1}='Original';
cond_str{2}='Null';
interestingConditions=[1,2];
%% Load Movies
rawMovFrames=960;
[stim,height,width,header_size] = get_raw_movie('/Volumes/Data/2015-01-29-2/Visual/18.rawMovie',rawMovFrames,1);
subtract_movies{3}=mean(stim,1);
subtract_movies{3}=mean(stim,1)*0+127.5;
movie=stim-repmat(subtract_movies{3},[rawMovFrames,1,1]);
movie=movie/255;

movStart=4*120*[0:6];
condMovies=cell(6,1);
for icond=1:2
    condMovies{icond}=movie(movStart(icond)+1:movStart(icond+1),:,:);
end

figure
for icond=1:2
    subplot(2,1,icond);
    qq=condMovies{icond}(120:4*120-120,:,:);
    hist(qq(:),20)
    title(sprintf('Condition %d',icond));
 end


figure('Color','w');
icnt=0;
for icondi=[2]
    icnt=icnt+1;
    subplot(1,1,icnt);
    aaa=condMovies{icondi}(120:120*4-120,:,:)-condMovies{1}(120:120*4-120,:,:);
    hist(aaa(:),50);
    title('Difference in pixel value compared to original movie');
end


%%


WN_datafile = '2015-01-29-2/data002_nps/data002_nps';
Null_datafile = '/Volumes/Analysis/2015-01-29-2/data005-from-data002_nps';
%WN_mov='/Volumes/Analysis/stimuli/white-noise-xml/RGB-10-2-0.48-11111.xml';
WN_datafile_full = '/Volumes/Analysis/2015-01-29-2/data002_nps';
imov=05;

datarun=load_data(WN_datafile)
datarun=load_params(datarun)



for cellTypeId=[1:10]; % 1 for On Parasols, 2 for Off parasols
InterestingCell_vis_id=[];
for icellType=cellTypeId  
    icellType 
    InterestingCell_vis_id=[InterestingCell_vis_id,datarun.cell_types{icellType}.cell_ids];
end 
cellTypeUsed=cellTypeId*ones(length(InterestingCell_vis_id),1);
NullCells=[2357,241,6211,4907,4160];
condDuration=3*6;

for ref_cell_number=1:length(InterestingCell_vis_id); %11
    close all
     [spkColl,spkCondColl,h]=plot_raster_script_2015_01_29_2(datarun,WN_datafile,WN_datafile_full,Null_datafile,InterestingCell_vis_id,imov,ref_cell_number,nConditions,condDuration,cond_str);
 
  plot_mosaic_2015_01_29_2(datarun,InterestingCell_vis_id,ref_cell_number,NullCells);
   InterestingCell_vis_id(ref_cell_number)
    if(~isdir(sprintf('/Volumes/Analysis/nishal/analyse_2015_01_29_2/CellType_%s',datarun.cell_types{cellTypeId}.name)))
    mkdir(sprintf('/Volumes/Analysis/nishal/analyse_2015_01_29_2/CellType_%s',datarun.cell_types{cellTypeId}.name));
    end
   s=hgexport('readstyle','ras_mos3');
    hgexport(h,sprintf('/Volumes/Analysis/nishal/analyse_2015_01_29_2/CellType_%s/CellID_%d.eps',datarun.cell_types{cellTypeId}.name,InterestingCell_vis_id(ref_cell_number)),s);
  
    %testsuite_prediction
  %[timeLogData,psthData] = psth_variability(spkCondColl,nConditions,condMovies,cond_str,InterestingCell_vis_id,imov,ref_cell_number,interestingConditions);
   % pause
end
end
