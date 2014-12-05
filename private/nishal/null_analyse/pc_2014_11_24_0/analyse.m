%%

addpath(genpath('../null_analyse/'));
addpath(genpath('../null_analyse/analyse_functions'));
startup_null_analyse_tenessee
%startup_null_analyse_bertha

%%
% Condition strings
nConditions=6;
condDuration=12;
cond_str=cell(6,1);
cond_str{1}='Original';
cond_str{2}='Null for On Parasol';
cond_str{3}='Original';
cond_str{4}='Null for few cells';
cond_str{5}='Original';
cond_str{6}='Null for Off Parasol';
interestingConditions=[1,2,4,6];
%% Load Movies

rawMovFrames=8640;
[stim,height,width,header_size] = get_raw_movie('/Volumes/Analysis/nishal/pc2014_11_24_0_data000/18.rawMovie',rawMovFrames,1);
subtract_movies{3}=mean(stim,1);
subtract_movies{3}=mean(stim,1)*0+127.5;
movie=stim-repmat(subtract_movies{3},[rawMovFrames,1,1]);
movie=movie/255;

movStart=12*120*[0:6];
condMovies=cell(6,1);
for icond=1:6
    condMovies{icond}=movie(movStart(icond)+1:movStart(icond+1),:,:);
end

figure
for icond=1:6
    subplot(3,2,icond);
    qq=condMovies{icond}(120:12*120-120,:,:);
    hist(qq(:),20)
    title(sprintf('Condition %d',icond));
end


figure('Color','w');
icnt=0;
for icondi=[2,4,6]
    icnt=icnt+1;
    subplot(3,1,icnt);
    aaa=condMovies{icondi}(240:1440-240,:,:)-condMovies{1}(240:1440-240,:,:);
    hist(aaa(:),50);
    title('Difference in pixel value compared to original movie');
end


%%



% On parasol
WN_datafile = 'nishal/2014-11-24-0/data000/data000';
Null_datafile = '/Volumes/Analysis/2014-11-24-0/data004';
WN_mov='/Volumes/Analysis/stimuli/white-noise-xml/RGB-10-2-0.48-11111.xml';
WN_datafile_full = '/Volumes/Analysis/nishal/2014-11-24-0/data000';
imov=04;

datarun=load_data(WN_datafile)
datarun=load_params(datarun)


cellTypeId=[1]; % 1 for On Parasols, 2 for Off parasols
InterestingCell_vis_id=[];
for icellType=cellTypeId
    icellType 
    InterestingCell_vis_id=[InterestingCell_vis_id,datarun.cell_types{icellType}.cell_ids];
end 
cellTypeUsed=cellTypeId*ones(length(InterestingCell_vis_id),1);

for ref_cell_number=[11,10,12,22,30]%:length(InterestingCell_vis_id); %11
    close all
    
    ref_cell_number
     [spkColl,spkCondColl]=plot_raster_script(datarun,WN_datafile,WN_datafile_full,Null_datafile,InterestingCell_vis_id,imov,ref_cell_number,nConditions,condDuration,cond_str);
      gen_clip=linear_output(WN_datafile,InterestingCell_vis_id,ref_cell_number,condMovies,'clip',nConditions,cond_str);
    gen_fit=linear_output(WN_datafile,InterestingCell_vis_id,ref_cell_number,condMovies,'fit',nConditions,cond_str);
    ref_cell_number
   
   plot_mosaic(datarun,InterestingCell_vis_id,ref_cell_number)
   
    %testsuite_prediction
  [timeLogData,psthData] = psth_variability(spkCondColl,nConditions,condMovies,cond_str,InterestingCell_vis_id,imov,ref_cell_number,interestingConditions);
    pause
end