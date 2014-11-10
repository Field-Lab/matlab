addpath(genpath('../null_analyse/'));
startup_null_analyse_tenessee


%% Analyse movies


rawMovFrames=5760;
[stim,height,width,header_size] = get_raw_movie('/Volumes/Data/2014-11-05-2/visual/18.rawMovie',rawMovFrames,1);
subtract_movies{3}=mean(stim,1);
subtract_movies{3}=mean(stim,1)*0+127.5;
movie=stim-repmat(subtract_movies{3},[rawMovFrames,1,1]);
movie=movie/255;

movStart=12*120*[0:4];
condMovies=cell(4,1);
for icond=1:4
condMovies{icond}=movie(movStart(icond)+1:movStart(icond+1),:,:);
end

figure
for icond=1:4
subplot(2,2,icond);
qq=condMovies{icond}(120:12*120-120,:,:);
hist(qq(:),20)
title(sprintf('Condition %d',icond));
end

pause
figure;
for itime=500:10:1440
    itime
    for icond=1:4
subplot(2,2,icond);
imagesc(squeeze(condMovies{icond}(itime,:,:)));
caxis([-0.5,0.5]);
colormap gray
axis image
colorbar

    end
pause
end

% Seems like the condition 3 (Off parasol) has lower contrast. This could
% be explained by the fact that having so manhy cells would have distorted
% pixel values by a lot. And as I tolerate only 2% of pixels being clipped
% off, the contrast had to decline!! 

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

 InterestingCell_vis_id=[3692,6061,1382,3782,6421,1627,2181,2631];


 for ref_cell_number=1:8
 plot_raster_script;
 pause
 end


