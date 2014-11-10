addpath(genpath('../null_analyse/'));
startup_null_analyse_tenessee


%%


rawMovFrames=5760;
[stim,height,width,header_size] = get_raw_movie('/Volumes/Data/2014-11-05-2/visual/18.rawMovie',rawMovFrames,1);
subtract_movies{3}=mean(stim,1);
subtract_movies{3}=mean(stim,1)*0+127.5;
movie=stim-repmat(subtract_movies{3},[rawMovFrames,1,1]);


%% Load cells and STAs

datafile = '2014-11-05-2/data009';
datarun=load_data(datafile)
datarun=load_params(datarun)


cellTypeId=[2]; % 1 for On Parasols, 2 for Off parasols
InterestingCell_vis_id=[];
for icellType=cellTypeId
    icellType
InterestingCell_vis_id=[InterestingCell_vis_id,datarun.cell_types{icellType}.cell_ids];
end

%%
 InterestingCell_vis_id=[3692,6061,1382,3782,6421,1627,2181];

 
 %%
% Off parasol
 datafile = '2014-11-05-2/data009';
datarun=load_data(datafile)
datarun=load_params(datarun)


cellTypeId=[2]; % 1 for On Parasols, 2 for Off parasols
InterestingCell_vis_id=[];
for icellType=cellTypeId
    icellType
InterestingCell_vis_id=[InterestingCell_vis_id,datarun.cell_types{icellType}.cell_ids];
end

 for ref_cell_number=1:10;
 plot_raster_script;
 pause
 end
%%
% On parasol
 datafile = '2014-11-05-2/data009';
datarun=load_data(datafile)
datarun=load_params(datarun)


cellTypeId=[1]; % 1 for On Parasols, 2 for Off parasols
InterestingCell_vis_id=[];
for icellType=cellTypeId
    icellType
InterestingCell_vis_id=[InterestingCell_vis_id,datarun.cell_types{icellType}.cell_ids];
end

 for ref_cell_number=1:10;
 plot_raster_script;
 pause
 end

 %%
% On parasol
 datafile = '2014-11-05-2/data009';
datarun=load_data(datafile)
datarun=load_params(datarun)

 InterestingCell_vis_id=[3692,6061,1382,3782,6421,1627,2181];


 for ref_cell_number=1:7;
 plot_raster_script;
 pause
 end


