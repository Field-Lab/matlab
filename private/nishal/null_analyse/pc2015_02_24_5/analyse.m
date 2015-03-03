addpath(genpath('../null_analyse/'));
addpath(genpath('../null_analyse/analyse_functions'));
%startup_null_analyse_tenessee
startup_null_analyse_bertha

%%
% Condition strings
nConditions=3;
condDuration=12;
cond_str=cell(3,1);
cond_str{1}='Original';
cond_str{2}='Null for cell group 1';
cond_str{3}='Null for cell group 2';
interestingConditions=[1,2,3];
%% Load Movies
rawMovFrames=4320/(8);
[stim,height,width,header_size] = get_raw_movie('/Volumes/Data/2015-02-24-5/Visual/pc2015_02_24_5_data001/19.rawMovie',rawMovFrames,1);
subtract_movies{3}=mean(stim,1);
subtract_movies{3}=mean(stim,1)*0+127.5;
movie=stim-repmat(subtract_movies{3},[rawMovFrames,1,1]);
movie=movie/255;

movStart=12*120*[0:3]/8;
condMovies=cell(3,1);
for icond=1:3
    condMovies{icond}=movie(movStart(icond)+1:movStart(icond+1),:,:);
end

figure
for icond=1:3
    subplot(3,2,icond);
    qq=condMovies{icond}(120/8:(12*120-120)/8,:,:);
    hist(qq(:),20)
    title(sprintf('Condition %d',icond));
end


figure('Color','w');
icnt=0;
for icondi=[1,2,3]
    icnt=icnt+1;
    subplot(3,1,icnt);
    aaa=condMovies{icondi}(240/8:(1440-240)/8,:,:)-condMovies{1}(240/8:(1440-240)/8,:,:);
    hist(aaa(:),50);
    title('Difference in pixel value compared to original movie');
end


%% data011


WN_datafile = '2015-02-24-5/streamed/data006/data006';
Null_datafile = '/Volumes/Analysis/2015-02-24-5/data011-from-data006_s_nps';
%WN_mov='/Volumes/Analysis/stimuli/white-noise-xml/RGB-10-2-0.48-11111.xml';
WN_datafile_full = '/Volumes/Analysis/2015-02-24-5/streamed/data006/data006';
neuronPath = [Null_datafile,sprintf('/data011-from-data006_s_nps.neurons')];
imov=14;

datarun=load_data(WN_datafile)
datarun=load_params(datarun)


cellTypeId=[13]; % 1 for On Parasols, 13 for ON large
InterestingCell_vis_id=[];
for icellType=cellTypeId
    icellType 
    InterestingCell_vis_id=[InterestingCell_vis_id,datarun.cell_types{icellType}.cell_ids];
end 

nConditions=2;
cellTypeUsed=cellTypeId*ones(length(InterestingCell_vis_id),1);
NullCells1=[6301,3181,2041,286];
NullCells2 = [3919];
condDuration=92;

for ref_cell_number=1:length(InterestingCell_vis_id); %11
    close all
     [spkColl,spkCondColl,h]=plot_raster_script_pc2015_02_24_5(datarun,WN_datafile,WN_datafile_full,Null_datafile,InterestingCell_vis_id,imov,ref_cell_number,nConditions,condDuration,cond_str,neuronPath);
 
  plot_mosaic_pc2015_02_24_5(datarun,InterestingCell_vis_id,ref_cell_number,NullCells1,NullCells2);
   InterestingCell_vis_id(ref_cell_number)
    if(~isdir(sprintf('/Volumes/Analysis/nishal/analyse_2015_02_24_5/data011/CellType_%s',datarun.cell_types{cellTypeId}.name)))
    mkdir(sprintf('/Volumes/Analysis/nishal/analyse_2015_02_24_5/data011/CellType_%s',datarun.cell_types{cellTypeId}.name));
    end
   s=hgexport('readstyle','ras_mos4');
   hgexport(h,sprintf('/Volumes/Analysis/nishal/analyse_2015_02_24_5/data011/CellType_%s/CellID_%d.eps',datarun.cell_types{cellTypeId}.name,InterestingCell_vis_id(ref_cell_number)),s);
  
    %testsuite_prediction
  %[timeLogData,psthData] = psth_variability(spkCondColl,nConditions,condMovies,cond_str,InterestingCell_vis_id,imov,ref_cell_number,interestingConditions);
   % pause
end

%% data010


WN_datafile = '2015-02-24-5/streamed/data006/data006';
Null_datafile = '/Volumes/Analysis/2015-02-24-5/data010-from-data006_s_nps';
%WN_mov='/Volumes/Analysis/stimuli/white-noise-xml/RGB-10-2-0.48-11111.xml';
WN_datafile_full = '/Volumes/Analysis/2015-02-24-5/streamed/data006/data006';
neuronPath = [Null_datafile,sprintf('/data010-from-data006_s_nps.neurons')];
imov=14;

datarun=load_data(WN_datafile)
datarun=load_params(datarun)


cellTypeId=[13]; % 1 for On Parasols, 13 for ON large
InterestingCell_vis_id=[];
for icellType=cellTypeId
    icellType 
    InterestingCell_vis_id=[InterestingCell_vis_id,datarun.cell_types{icellType}.cell_ids];
end 

nConditions=3;
cellTypeUsed=cellTypeId*ones(length(InterestingCell_vis_id),1);
NullCells1=[6301,3181,2041,286];
NullCells2 = [3919];
condDuration=12;

for ref_cell_number=1:length(InterestingCell_vis_id); %11
    close all
     [spkColl,spkCondColl,h]=plot_raster_script_pc2015_02_24_5(datarun,WN_datafile,WN_datafile_full,Null_datafile,InterestingCell_vis_id,imov,ref_cell_number,nConditions,condDuration,cond_str,neuronPath);
 
  plot_mosaic_pc2015_02_24_5(datarun,InterestingCell_vis_id,ref_cell_number,NullCells1,NullCells2);
   InterestingCell_vis_id(ref_cell_number)
    if(~isdir(sprintf('/Volumes/Analysis/nishal/analyse_2015_02_24_5/data010/CellType_%s',datarun.cell_types{cellTypeId}.name)))
    mkdir(sprintf('/Volumes/Analysis/nishal/analyse_2015_02_24_5/data010/CellType_%s',datarun.cell_types{cellTypeId}.name));
    end
   s=hgexport('readstyle','ras_mos4');
   hgexport(h,sprintf('/Volumes/Analysis/nishal/analyse_2015_02_24_5/data010/CellType_%s/CellID_%d.eps',datarun.cell_types{cellTypeId}.name,InterestingCell_vis_id(ref_cell_number)),s);
  
    %testsuite_prediction
  %[timeLogData,psthData] = psth_variability(spkCondColl,nConditions,condMovies,cond_str,InterestingCell_vis_id,imov,ref_cell_number,interestingConditions);
   % pause
end




%% data004


WN_datafile = '2015-02-24-5/streamed/data001/data001';
Null_datafile = '/Volumes/Analysis/2015-02-24-5/data004-from-data001_s_nps';
%WN_mov='/Volumes/Analysis/stimuli/white-noise-xml/RGB-10-2-0.48-11111.xml';
WN_datafile_full = '/Volumes/Analysis/2015-02-24-5/streamed/data001/data001';
neuronPath = [Null_datafile,sprintf('/data004-from-data001_s_nps.neurons')];
imov=14;

datarun=load_data(WN_datafile)
datarun=load_params(datarun)


cellTypeId=[10]; % 1 for On Parasols, 8 for OFF, 10, for ON nc2
InterestingCell_vis_id=[];
for icellType=cellTypeId
    icellType 
    InterestingCell_vis_id=[InterestingCell_vis_id,datarun.cell_types{icellType}.cell_ids];
end 

nConditions=3;
cellTypeUsed=cellTypeId*ones(length(InterestingCell_vis_id),1);
NullCells1=[6121,3181,286,3722];
NullCells2 = [7518,1306,3723];
condDuration=12;

for ref_cell_number=1:length(InterestingCell_vis_id); %11
    close all
     [spkColl,spkCondColl,h]=plot_raster_script_pc2015_02_24_5(datarun,WN_datafile,WN_datafile_full,Null_datafile,InterestingCell_vis_id,imov,ref_cell_number,nConditions,condDuration,cond_str,neuronPath);
 
  plot_mosaic_pc2015_02_24_5(datarun,InterestingCell_vis_id,ref_cell_number,NullCells1,NullCells2);
   InterestingCell_vis_id(ref_cell_number)
    if(~isdir(sprintf('/Volumes/Analysis/nishal/analyse_2015_02_24_5/data004/CellType_%s',datarun.cell_types{cellTypeId}.name)))
    mkdir(sprintf('/Volumes/Analysis/nishal/analyse_2015_02_24_5/data004/CellType_%s',datarun.cell_types{cellTypeId}.name));
    end
   s=hgexport('readstyle','ras_mos4');
   hgexport(h,sprintf('/Volumes/Analysis/nishal/analyse_2015_02_24_5/data004/CellType_%s/CellID_%d.eps',datarun.cell_types{cellTypeId}.name,InterestingCell_vis_id(ref_cell_number)),s);
  
    %testsuite_prediction
  %[timeLogData,psthData] = psth_variability(spkCondColl,nConditions,condMovies,cond_str,InterestingCell_vis_id,imov,ref_cell_number,interestingConditions);
   % pause
end


%% data005


WN_datafile = '2015-02-24-5/streamed/data001/data001';
Null_datafile = '/Volumes/Analysis/2015-02-24-5/data005-from-data001_s_nps';
%WN_mov='/Volumes/Analysis/stimuli/white-noise-xml/RGB-10-2-0.48-11111.xml';
WN_datafile_full = '/Volumes/Analysis/2015-02-24-5/streamed/data001/data001';
neuronPath = [Null_datafile,sprintf('/data005-from-data001_s_nps.neurons')];
imov=14;

datarun=load_data(WN_datafile)
datarun=load_params(datarun)


cellTypeId=[10]; % 1 for On Parasols, 8 for OFF, 10 for ON nc2
InterestingCell_vis_id=[];
for icellType=cellTypeId
    icellType 
    InterestingCell_vis_id=[InterestingCell_vis_id,datarun.cell_types{icellType}.cell_ids];
end 

nConditions=2;
cellTypeUsed=cellTypeId*ones(length(InterestingCell_vis_id),1);
NullCells1=[6121,3181,286,3722];
NullCells2 = [7518,1306,3723];
condDuration=92;

for ref_cell_number=1:length(InterestingCell_vis_id); %11
    close all
     [spkColl,spkCondColl,h]=plot_raster_script_pc2015_02_24_5(datarun,WN_datafile,WN_datafile_full,Null_datafile,InterestingCell_vis_id,imov,ref_cell_number,nConditions,condDuration,cond_str,neuronPath);
 
  plot_mosaic_pc2015_02_24_5(datarun,InterestingCell_vis_id,ref_cell_number,NullCells1,NullCells2);
   InterestingCell_vis_id(ref_cell_number)
    if(~isdir(sprintf('/Volumes/Analysis/nishal/analyse_2015_02_24_5/data005/CellType_%s',datarun.cell_types{cellTypeId}.name)))
    mkdir(sprintf('/Volumes/Analysis/nishal/analyse_2015_02_24_5/data005/CellType_%s',datarun.cell_types{cellTypeId}.name));
    end
   s=hgexport('readstyle','ras_mos4');
   hgexport(h,sprintf('/Volumes/Analysis/nishal/analyse_2015_02_24_5/data005/CellType_%s/CellID_%d.eps',datarun.cell_types{cellTypeId}.name,InterestingCell_vis_id(ref_cell_number)),s);
  
    %testsuite_prediction
  %[timeLogData,psthData] = psth_variability(spkCondColl,nConditions,condMovies,cond_str,InterestingCell_vis_id,imov,ref_cell_number,interestingConditions);
   % pause
end



%% data020


WN_datafile = '2015-02-24-5/streamed/data018/data018';
Null_datafile = '/Volumes/Analysis/2015-02-24-5/data020-from-data018_s_nps';
%WN_mov='/Volumes/Analysis/stimuli/white-noise-xml/RGB-10-2-0.48-11111.xml';
WN_datafile_full = '/Volumes/Analysis/2015-02-24-5/streamed/data018/data018';
neuronPath = [Null_datafile,sprintf('/data020-from-data018_s_nps.neurons')];
imov=14;

datarun=load_data(WN_datafile)
datarun=load_params(datarun)


cellTypeId=[12]; % 17 for On Parasols, 8 for large,12 for OFF
InterestingCell_vis_id=[];
for icellType=cellTypeId
    icellType 
    InterestingCell_vis_id=[InterestingCell_vis_id,datarun.cell_types{icellType}.cell_ids];
end 

nConditions=3;
cellTypeUsed=cellTypeId*ones(length(InterestingCell_vis_id),1);
NullCells1='ON';
NullCells2 =[2432,16,4996];
condDuration=12;

for ref_cell_number=1:length(InterestingCell_vis_id); %11
    close all
     [spkColl,spkCondColl,h]=plot_raster_script_pc2015_02_24_5(datarun,WN_datafile,WN_datafile_full,Null_datafile,InterestingCell_vis_id,imov,ref_cell_number,nConditions,condDuration,cond_str,neuronPath);
 
  plot_mosaic_pc2015_02_24_5(datarun,InterestingCell_vis_id,ref_cell_number,NullCells1,NullCells2);
   InterestingCell_vis_id(ref_cell_number)
    if(~isdir(sprintf('/Volumes/Analysis/nishal/analyse_2015_02_24_5/data020/CellType_%s',datarun.cell_types{cellTypeId}.name)))
    mkdir(sprintf('/Volumes/Analysis/nishal/analyse_2015_02_24_5/data020/CellType_%s',datarun.cell_types{cellTypeId}.name));
    end
   s=hgexport('readstyle','ras_mos4');
   hgexport(h,sprintf('/Volumes/Analysis/nishal/analyse_2015_02_24_5/data020/CellType_%s/CellID_%d.eps',datarun.cell_types{cellTypeId}.name,InterestingCell_vis_id(ref_cell_number)),s);
  
    %testsuite_prediction
  %[timeLogData,psthData] = psth_variability(spkCondColl,nConditions,condMovies,cond_str,InterestingCell_vis_id,imov,ref_cell_number,interestingConditions);
   % pause
end


%% data021


WN_datafile = '2015-02-24-5/streamed/data018/data018';
Null_datafile = '/Volumes/Analysis/2015-02-24-5/data021-from-data018_s_nps';
%WN_mov='/Volumes/Analysis/stimuli/white-noise-xml/RGB-10-2-0.48-11111.xml';
WN_datafile_full = '/Volumes/Analysis/2015-02-24-5/streamed/data018/data018';
neuronPath = [Null_datafile,sprintf('/data021-from-data018_s_nps.neurons')];
imov=14;

datarun=load_data(WN_datafile)
datarun=load_params(datarun)


cellTypeId=[12]; % 17 for On Parasols, 8 for large,12 for OFF
InterestingCell_vis_id=[];
for icellType=cellTypeId
    icellType 
    InterestingCell_vis_id=[InterestingCell_vis_id,datarun.cell_types{icellType}.cell_ids];
end 

nConditions=3;
cellTypeUsed=cellTypeId*ones(length(InterestingCell_vis_id),1);
NullCells1='ON';
NullCells2 =[2432,16,4996];
condDuration=12;

for ref_cell_number=1:length(InterestingCell_vis_id); %11
    close all
     [spkColl,spkCondColl,h]=plot_raster_script_pc2015_02_24_5(datarun,WN_datafile,WN_datafile_full,Null_datafile,InterestingCell_vis_id,imov,ref_cell_number,nConditions,condDuration,cond_str,neuronPath);
 
  plot_mosaic_pc2015_02_24_5(datarun,InterestingCell_vis_id,ref_cell_number,NullCells1,NullCells2);
   InterestingCell_vis_id(ref_cell_number)
    if(~isdir(sprintf('/Volumes/Analysis/nishal/analyse_2015_02_24_5/data021/CellType_%s',datarun.cell_types{cellTypeId}.name)))
    mkdir(sprintf('/Volumes/Analysis/nishal/analyse_2015_02_24_5/data021/CellType_%s',datarun.cell_types{cellTypeId}.name));
    end
   s=hgexport('readstyle','ras_mos4');
   hgexport(h,sprintf('/Volumes/Analysis/nishal/analyse_2015_02_24_5/data021/CellType_%s/CellID_%d.eps',datarun.cell_types{cellTypeId}.name,InterestingCell_vis_id(ref_cell_number)),s);
  
    %testsuite_prediction
  %[timeLogData,psthData] = psth_variability(spkCondColl,nConditions,condMovies,cond_str,InterestingCell_vis_id,imov,ref_cell_number,interestingConditions);
   % pause
end


