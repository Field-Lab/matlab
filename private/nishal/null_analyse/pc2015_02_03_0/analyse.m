addpath(genpath('../null_analyse/'));
addpath(genpath('../null_analyse/analyse_functions'));
%startup_null_analyse_tenessee
startup_null_analyse_bertha

%%
% Condition strings
nConditions=5;
condDuration=12;
cond_str=cell(nConditions,1);
cond_str{1}='Original';
cond_str{2}='Spatial Null for Group 1';
cond_str{3}='Spatio-Temporal Null for Group 1';
cond_str{4}='Spatial Null for Group 2';
cond_str{5}='Spatio-Temporal Null for Group 2';

interestingConditions=[1,2,3,4,5];
%% Load Movies
rawMovFrames=7200;
[stim,height,width,header_size] = get_raw_movie('/Volumes/Data/2015-02-03-0/Visual/18.rawMovie',rawMovFrames,1);

subtract_movies=mean(stim,1)*0+127.5;
movie=stim-repmat(subtract_movies,[rawMovFrames,1,1]);
movie=movie/255;

movStart=12*120*[0:6];
condMovies=cell(6,1);
for icond=1:nConditions
    condMovies{icond}=movie(movStart(icond)+1:movStart(icond+1),:,:);
end

figure
for icond=1:nConditions
    subplot(nConditions,1,icond);
    qq=condMovies{icond}(120:4*120-120,:,:);
    hist(qq(:),20)
    title(sprintf('Condition %d',icond));
 end


figure('Color','w');
icnt=0;
for icondi=[2,3,4,5]
    icnt=icnt+1;
    subplot(4,1,icnt);
    aaa=condMovies{icondi}(120:120*4-120,:,:)-condMovies{1}(120:120*4-120,:,:);
    hist(aaa(:),50);
    title(sprintf('Difference in pixel value compared to original movie: Cond %d',icondi));
end


%%


WN_datafile = '2015-02-03-0/data000/data000';
Null_datafile = '/Volumes/Analysis/2015-02-03-0/data001-from-data000';
%WN_mov='/Volumes/Analysis/stimuli/white-noise-xml/RGB-10-2-0.48-11111.xml';
WN_datafile_full = '/Volumes/Analysis/2015-02-03-0/data000';
imov=01;

datarun=load_data(WN_datafile)
datarun=load_params(datarun)



for cellTypeId=[1:18]; % 1 for On Parasols, 2 for Off parasols
    try
InterestingCell_vis_id=[];
for icellType=cellTypeId  
    icellType 
    InterestingCell_vis_id=[InterestingCell_vis_id,datarun.cell_types{icellType}.cell_ids];
end 
cellTypeUsed=cellTypeId*ones(length(InterestingCell_vis_id),1);
NullCells1=[6181,7022,6841,4966,4876,5416,818,976,1187,1205,1488,1853,4071,4727,4831,5131,6182];
NullCells2 = [4876,4966,5416,6181,6841,7022];
condDuration=12;

for ref_cell_number=1:length(InterestingCell_vis_id); %11
    close all
     [spkColl,spkCondColl,h]=plot_raster_script_pc2015_02_03_0(datarun,WN_datafile,WN_datafile_full,Null_datafile,InterestingCell_vis_id,imov,ref_cell_number,nConditions,condDuration,cond_str);
 
  plot_mosaic_pc2015_02_03_0(datarun,InterestingCell_vis_id,ref_cell_number,NullCells1,NullCells2);
   InterestingCell_vis_id(ref_cell_number)
    if(~isdir(sprintf('/Volumes/Analysis/nishal/analyse_2015_02_03_0/CellType_%s',datarun.cell_types{cellTypeId}.name)))
    mkdir(sprintf('/Volumes/Analysis/nishal/analyse_2015_02_03_0/CellType_%s',datarun.cell_types{cellTypeId}.name));
    end
   s=hgexport('readstyle','ras_mos4');
    hgexport(h,sprintf('/Volumes/Analysis/nishal/analyse_2015_02_03_0/CellType_%s/CellID_%d.eps',datarun.cell_types{cellTypeId}.name,InterestingCell_vis_id(ref_cell_number)),s);
  
    %testsuite_prediction
  %[timeLogData,psthData] = psth_variability(spkCondColl,nConditions,condMovies,cond_str,InterestingCell_vis_id,imov,ref_cell_number,interestingConditions);
   % pause
end
    catch
    end
end
