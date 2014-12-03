addpath(genpath('../null_analyse/'));
startup_null_analyse_tenessee
%startup_null_analyse_bertha

%% Analyse movies


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



pause
figure;
for itime=500:10:1440
    itime
    for icond=1:6
        subplot(3,2,icond);
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

%%
% Condition strings
cond_str=cell(6,1);
cond_str{1}='Original';
cond_str{2}='Null for On Parasol';
cond_str{3}='Original';
cond_str{4}='Null for Off Parasol';
cond_str{5}='Original';
cond_str{6}='Null for few cells';
%% Load cells and STAs

WN_datafile = 'nishal/2014-11-24-0/data000/data000';

imov=04;
datarun=load_data(WN_datafile)
datarun=load_params(datarun)


cellTypeId=[1]; % 1 for On Parasols, 2 for Off parasols
InterestingCell_vis_id=[];
for icellType=cellTypeId
    icellType
    InterestingCell_vis_id=[InterestingCell_vis_id,datarun.cell_types{icellType}.cell_ids];
end

%%

% On parasol
WN_datafile = 'nishal/2014-11-24-0/data000/data000';
Null_datafile = '/Volumes/Analysis/2014-11-24-0/data004';
datarun=load_data(WN_datafile)
datarun=load_params(datarun)


cellTypeId=[8]; % 1 for On Parasols, 2 for Off parasols
InterestingCell_vis_id=[];
for icellType=cellTypeId
    icellType 
    InterestingCell_vis_id=[InterestingCell_vis_id,datarun.cell_types{icellType}.cell_ids];
end 
cellTypeUsed=cellTypeId*ones(length(InterestingCell_vis_id),1);
WN_datafile_full = '/Volumes/Analysis/nishal/2014-11-24-0/data000';
for ref_cell_number=1:length(InterestingCell_vis_id); %11
    close all
    
    ref_cell_number
    plot_raster_script;
    %plot_mosaic
    %testsuite_prediction
    %psth_variability
    pause
end


