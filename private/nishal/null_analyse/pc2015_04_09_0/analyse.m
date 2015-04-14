addpath(genpath('../null_analyse/'));
addpath(genpath('../null_analyse/analyse_functions'));
%startup_null_analyse_tenessee
startup_null_analyse_bertha

%%
% Condition strings
nConditions=5;
condDuration=1200/1;
cond_str=cell(5,1);
cond_str{1}='Original';
cond_str{2}='Null Spatial ';
cond_str{3}='WN_NULL Spatial';
cond_str{4}='Null Spatio-temporal ';
cond_str{5}='WN+Null Spatio-temporal';

interestingConditions=[1,2,3,4,5,6];
%% Load Movies

rawMovFrames=1200/(1);
figure;
icnt=0;
% make pixel histogram
for imov=[1,2,4,6,8]
    [stim,height,width,header_size] = get_raw_movie(sprintf('/Volumes/Data/2015-04-09-0/Visual/Null/pc2015_04_09_0_data000/%d.rawMovie',imov),rawMovFrames,1);
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

% make movies
interval=1;
condMov=cell(nConditions,1);
rawMovFrames=1200/(1);
icnt=0;
% make pixel histogram
for imov=[1,2,4,6,8]
    [stim,height,width,header_size] = get_raw_movie(sprintf('/Volumes/Data/2015-04-09-0/Visual/Null/pc2015_04_09_0_data000/%d.rawMovie',imov),rawMovFrames,1);
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
rawMovFrames=1200/(1);
figure;
icnt=0;
cMap = cell(6,1);
h=figure('Color','w');
for imov=[1,2,4,6,8]
    [stim,height,width,header_size] = get_raw_movie(sprintf('/Volumes/Data/2015-04-09-0/Visual/Null/pc2015_04_09_0_data000/%d.rawMovie',imov),rawMovFrames,1);
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




%% data006 from data000


WN_datafile = '2015-04-09-0/stix8-norefit/data000-from-data000_data003_data004_data005_data008_data009_data016_data017_data006/data000-from-data000_data003_data004_data005_data008_data009_data016_data017_data006';
Null_datafile = '/Volumes/Analysis/2015-04-09-0/stix8-norefit/data006-from-data000_data003_data004_data005_data008_data009_data016_data017_data006';
%WN_mov='/Volumes/Analysis/stimuli/white-noise-xml/RGB-10-2-0.48-11111.xml';
neuronPath = [Null_datafile,sprintf('/data006-from-data000_data003_data004_data005_data008_data009_data016_data017_data006.neurons')];
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
NullCells1=[4775,6757,7462,1685,2722,5948,1459]; % cell : 6826 ?



condDuration=10;
nConditions=5;
for ref_cell_number=1:length(InterestingCell_vis_id); %11
    close all
    cellID=InterestingCell_vis_id(ref_cell_number);
    % [spkColl,spkCondColl,h]=plot_raster_script_pc2015_02_24_2(datarun,WN_datafile,WN_datafile_full,Null_datafile,InterestingCell_vis_id,imov,ref_cell_number,nConditions,condDuration,cond_str,neuronPath);
    [spkColl,spkCondColl,h]=plot_raster_script_pc2015_04_09_0_light(cellID,nConditions,condDuration,cond_str,neuronPath);
    
    plot_mosaic_pc2015_04_09_0(datarun,InterestingCell_vis_id,ref_cell_number,NullCells1);
    InterestingCell_vis_id(ref_cell_number)
    if(~isdir(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_04_09_0/data006/CellType_%s',datarun.cell_types{cellTypeId}.name)))
        mkdir(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_04_09_0/data006/CellType_%s',datarun.cell_types{cellTypeId}.name));
    end
    s=hgexport('readstyle','ras_mos4');
    hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_04_09_0/data006/CellType_%s/CellID_%d.eps',datarun.cell_types{cellTypeId}.name,InterestingCell_vis_id(ref_cell_number)),s);
    
    %testsuite_prediction
    %[timeLogData,psthData] = psth_variability(spkCondColl,nConditions,condMovies,cond_str,InterestingCell_vis_id,imov,ref_cell_number,interestingConditions);
    % pause
end


%% data011 from data007


WN_datafile = '2015-04-09-0/stix16-norefit/data007-from-data010_data013_data014_data015_data011_data007_data012/data007-from-data010_data013_data014_data015_data011_data007_data012';
Null_datafile = '/Volumes/Analysis/2015-04-09-0/stix16-norefit/data011-from-data010_data013_data014_data015_data011_data007_data012';
%WN_mov='/Volumes/Analysis/stimuli/white-noise-xml/RGB-10-2-0.48-11111.xml';
neuronPath = [Null_datafile,sprintf('/data011-from-data010_data013_data014_data015_data011_data007_data012.neurons')];
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
NullCells1=[]; % cell : 6826 ?



condDuration=10;
nConditions=5;
for ref_cell_number=1:length(InterestingCell_vis_id); %11
    close all
    cellID=InterestingCell_vis_id(ref_cell_number);
    % [spkColl,spkCondColl,h]=plot_raster_script_pc2015_02_24_2(datarun,WN_datafile,WN_datafile_full,Null_datafile,InterestingCell_vis_id,imov,ref_cell_number,nConditions,condDuration,cond_str,neuronPath);
    [spkColl,spkCondColl,h]=plot_raster_script_pc2015_04_09_0_light(cellID,nConditions,condDuration,cond_str,neuronPath);
    
    plot_mosaic_pc2015_04_09_0(datarun,InterestingCell_vis_id,ref_cell_number,NullCells1);
    InterestingCell_vis_id(ref_cell_number)
    if(~isdir(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_04_09_0/data011/CellType_%s',datarun.cell_types{cellTypeId}.name)))
        mkdir(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_04_09_0/data011/CellType_%s',datarun.cell_types{cellTypeId}.name));
    end
    s=hgexport('readstyle','ras_mos4');
    hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_04_09_0/data011/CellType_%s/CellID_%d.eps',datarun.cell_types{cellTypeId}.name,InterestingCell_vis_id(ref_cell_number)),s);
    
    %testsuite_prediction
    %[timeLogData,psthData] = psth_variability(spkCondColl,nConditions,condMovies,cond_str,InterestingCell_vis_id,imov,ref_cell_number,interestingConditions);
    % pause
end
