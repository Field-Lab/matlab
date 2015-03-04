addpath(genpath('../null_analyse/'));
addpath(genpath('../null_analyse/analyse_functions'));
startup_null_analyse_tenessee
%startup_null_analyse_bertha

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
[stim,height,width,header_size] = get_raw_movie('/Volumes/Data/2015-02-24-2/Visual/pc_2015_02_24_2_data010/19.rawMovie',rawMovFrames,1);
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

% Contrast map
cMap = cell(3,1);
h=figure('Color','w');
for icond=1:3
    cMap{icond}=contrastMap(condMovies{icond});
    subplot(2,2,icond);
    imagesc(cMap{icond});
    caxis([3,6]);
    colorbar
    title(sprintf('cMap: %s',cond_str{icond}));
end

   s=hgexport('readstyle','cMap');
   hgexport(h,sprintf('/Volumes/Analysis/nishal/analyse_2015_02_24_2/data014/cMap.eps'),s);
  
%%


WN_datafile = '2015-02-24-2/streamed/data006/data006';
Null_datafile = '/Volumes/Analysis/2015-02-24-2/data012-from-data006_s_nps';
%WN_mov='/Volumes/Analysis/stimuli/white-noise-xml/RGB-10-2-0.48-11111.xml';
WN_datafile_full = '/Volumes/Analysis/2015-02-24-2/streamed/data006/data006';
neuronPath = [Null_datafile,sprintf('/data012-from-data006_s_nps.neurons')];
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
NullCells1=[5269,3706,2417,5192,7651,5897];
NullCells2 = [6211,1562,2971,5851,1383,5779];
condDuration=12;

for ref_cell_number=1:length(InterestingCell_vis_id); %11
    close all
     [spkColl,spkCondColl,h]=plot_raster_script_pc2015_02_24_2(datarun,WN_datafile,WN_datafile_full,Null_datafile,InterestingCell_vis_id,imov,ref_cell_number,nConditions,condDuration,cond_str,neuronPath);
 
  plot_mosaic_pc2015_02_24_2(datarun,InterestingCell_vis_id,ref_cell_number,NullCells1,NullCells2);
   InterestingCell_vis_id(ref_cell_number)
    if(~isdir(sprintf('/Volumes/Analysis/nishal/analyse_2015_02_24_2/data012/CellType_%s',datarun.cell_types{cellTypeId}.name)))
    mkdir(sprintf('/Volumes/Analysis/nishal/analyse_2015_02_24_2/data012/CellType_%s',datarun.cell_types{cellTypeId}.name));
    end
   s=hgexport('readstyle','ras_mos4');
   hgexport(h,sprintf('/Volumes/Analysis/nishal/analyse_2015_02_24_2/data012/CellType_%s/CellID_%d.eps',datarun.cell_types{cellTypeId}.name,InterestingCell_vis_id(ref_cell_number)),s);
  
    %testsuite_prediction
  %[timeLogData,psthData] = psth_variability(spkCondColl,nConditions,condMovies,cond_str,InterestingCell_vis_id,imov,ref_cell_number,interestingConditions);
   % pause
end

%
WN_datafile = '2015-02-24-2/streamed/data016/data016';
Null_datafile = '/Volumes/Analysis/2015-02-24-2/data019-from-data016_s_nps';
%WN_mov='/Volumes/Analysis/stimuli/white-noise-xml/RGB-10-2-0.48-11111.xml';
WN_datafile_full = '/Volumes/Analysis/2015-02-24-2/streamed/data016/data016';
neuronPath = [Null_datafile,sprintf('/data019-from-data016_s_nps.neurons')];
imov=19;

datarun=load_data(WN_datafile)
datarun=load_params(datarun)


cellTypeId=[1]; % 1 for On Parasols, 2 for Off parasols
InterestingCell_vis_id=[];
for icellType=cellTypeId
    icellType 
    InterestingCell_vis_id=[InterestingCell_vis_id,datarun.cell_types{icellType}.cell_ids];
end 

cellTypeUsed=cellTypeId*ones(length(InterestingCell_vis_id),1);
NullCells1=[2793,1832,5131,6948];
NullCells2 = [5042,1531,2911,7246];
condDuration=12;

for ref_cell_number=1:length(InterestingCell_vis_id); %11
    close all
     [spkColl,spkCondColl,h]=plot_raster_script_pc2015_02_24_2(datarun,WN_datafile,WN_datafile_full,Null_datafile,InterestingCell_vis_id,imov,ref_cell_number,nConditions,condDuration,cond_str,neuronPath);
 
  plot_mosaic_pc2015_02_24_2(datarun,InterestingCell_vis_id,ref_cell_number,NullCells1,NullCells2);
   InterestingCell_vis_id(ref_cell_number)
    if(~isdir(sprintf('/Volumes/Analysis/nishal/analyse_2015_02_24_2/data016/CellType_%s',datarun.cell_types{cellTypeId}.name)))
    mkdir(sprintf('/Volumes/Analysis/nishal/analyse_2015_02_24_2/data016/CellType_%s',datarun.cell_types{cellTypeId}.name));
    end
   s=hgexport('readstyle','ras_mos4');
   hgexport(h,sprintf('/Volumes/Analysis/nishal/analyse_2015_02_24_2/data016/CellType_%s/CellID_%d.eps',datarun.cell_types{cellTypeId}.name,InterestingCell_vis_id(ref_cell_number)),s);
  
    %testsuite_prediction
  %[timeLogData,psthData] = psth_variability(spkCondColl,nConditions,condMovies,cond_str,InterestingCell_vis_id,imov,ref_cell_number,interestingConditions);
   % pause
end


%% Joint analysis - from 016

WN_datafile = '2015-02-24-2/streamed/data016/data016';
datarun1=load_data(WN_datafile);
datarun1=load_params(datarun1);
datarun1=load_ei(datarun1,datarun1.cell_ids);

WN_datafile_concat = '2015-02-24-2/dNull-norefit_nps/data016-from-data012_data014_data019_data006_data010_data016/data016-from-data012_data014_data019_data006_data010_data016';
datarun2=load_data(WN_datafile_concat);
datarun2=load_params(datarun2);
datarun2.names.rrs_ei_path='/Volumes/Analysis/2015-02-24-2/dNull-norefit_nps/dNull-norefit_nps.ei';
datarun2=load_ei(datarun2,datarun2.cell_ids);

[cell_list_map, failed_cells] =map_ei(datarun1,datarun2);

NullCells1_1=[2793,1832,5131,6948];
NullCells2_1 = [5042,1531,2911,7246]; 

NullCells1_2=[];
for iList=1:length(NullCells1_1)
x=cell_list_map(datarun1.cell_ids==NullCells1_1(iList));
NullCells1_2(iList)=x{1};
end

NullCells2_2=[];
for iList=1:length(NullCells2_1)
x=cell_list_map(datarun1.cell_ids==NullCells2_1(iList));
NullCells2_2(iList)=x{1};
end

%
cellListOrig=NullCells2_1;
cellList=NullCells2_2;
for icell=1:length(cellList);
    InterestingCell_vis_id=cellList(icell);


% Print data019
Null_datafile = '/Volumes/Analysis/2015-02-24-2/dNull-norefit_nps/data019-from-data012_data014_data019_data006_data010_data016';
neuronPath = [Null_datafile,sprintf('/data019-from-data012_data014_data019_data006_data010_data016.neurons')];


ref_cell_number=1;
nConditions=3;
condDuration=12;

[spkColl,spkCondColl019,h]=plot_raster_script_pc2015_02_24_2_light(InterestingCell_vis_id(ref_cell_number),nConditions,condDuration,cond_str,neuronPath)

% Print data014
Null_datafile = '/Volumes/Analysis/2015-02-24-2/dNull-norefit_nps/data014-from-data012_data014_data019_data006_data010_data016';
neuronPath = [Null_datafile,sprintf('/data014-from-data012_data014_data019_data006_data010_data016.neurons')];

imov=1;
ref_cell_number=1;
nConditions=3;
condDuration=12;


[spkColl,spkCondColl014,h]=plot_raster_script_pc2015_02_24_2_light(InterestingCell_vis_id(ref_cell_number),nConditions,condDuration,cond_str,neuronPath)

% Print data012
Null_datafile = '/Volumes/Analysis/2015-02-24-2/dNull-norefit_nps/data012-from-data012_data014_data019_data006_data010_data016';
neuronPath = [Null_datafile,sprintf('/data012-from-data012_data014_data019_data006_data010_data016.neurons')];

imov=1;
ref_cell_number=1;
nConditions=3;
condDuration=12;

[spkColl,spkCondColl012,h]=plot_raster_script_pc2015_02_24_2_light(InterestingCell_vis_id(ref_cell_number),nConditions,condDuration,cond_str,neuronPath)

col='rkmrkmrkm';
databaseShift = 3*30;
nTrials1=30;
h=figure('Color','w');
cnt=0;

for icond=1:3
    xPoints=spkCondColl019(icond).xPoints;
    yPoints=spkCondColl019(icond).yPoints;
plot(xPoints*1/20000, yPoints+cnt,'r');
cnt=cnt+max(yPoints);
hold on

    xPoints=spkCondColl014(icond).xPoints;
    yPoints=spkCondColl014(icond).yPoints;
plot(xPoints*1/20000, yPoints+cnt,'k');
cnt=cnt+max(yPoints);
hold on

    xPoints=spkCondColl012(icond).xPoints;
    yPoints=spkCondColl012(icond).yPoints;
plot(xPoints*1/20000, yPoints+cnt,'m');
cnt=cnt+max(yPoints);
hold on
end

xlim([0,condDuration]);
ylim([0,cnt]);
title('data012,014,019 - order: Null2,Null1,Original');

    if(~isdir(sprintf('/Volumes/Analysis/nishal/analyse_2015_02_24_2/data016/NullCells2')))
    mkdir(sprintf('/Volumes/Analysis/nishal/analyse_2015_02_24_2/data016/NullCells2'));
    end
   s=hgexport('readstyle','ras5');
   hgexport(h,sprintf('/Volumes/Analysis/nishal/analyse_2015_02_24_2/data016/NullCells2/CellID_%d.eps',cellListOrig(icell)),s);
  


end


%% Joint analysis - from 006

WN_datafile = '2015-02-24-2/streamed/data006/data006';
datarun1=load_data(WN_datafile);
datarun1=load_params(datarun1);
datarun1=load_ei(datarun1,datarun1.cell_ids);

WN_datafile_concat = '2015-02-24-2/dNull-norefit_nps/data006-from-data012_data014_data019_data006_data010_data016/data006-from-data012_data014_data019_data006_data010_data016';
datarun2=load_data(WN_datafile_concat);
datarun2=load_params(datarun2);
datarun2.names.rrs_ei_path='/Volumes/Analysis/2015-02-24-2/dNull-norefit_nps/dNull-norefit_nps.ei';
datarun2=load_ei(datarun2,datarun2.cell_ids);

[cell_list_map, failed_cells] =map_ei(datarun1,datarun2);

NullCells1=[5269,3706,2417,5192,7651,5897];
NullCells2 = [6211,1562,2971,5851,1383,5779];

NullCells1_2=[];
NullCells1_1=[];
icnt=1;
for iList=1:length(NullCells1)
x=cell_list_map(datarun1.cell_ids==NullCells1(iList));
try
NullCells1_2(icnt)=x{1};
NullCells1_1(icnt)=NullCells1(icnt);
icnt=icnt+1;
catch
display('No match');
end

end

NullCells2_2=[];
NullCells2_1=[];
icnt=1;
for iList=1:length(NullCells2)
x=cell_list_map(datarun1.cell_ids==NullCells2(iList));
try
NullCells2_2(icnt)=x{1};
NullCells2_1(icnt)=NullCells2(icnt);
icnt=icnt+1;
catch
display('No match');    
end
end

%
cellListOrig=NullCells1_1;
cellList=NullCells1_2;
for icell=1:length(cellList);
    InterestingCell_vis_id=cellList(icell);


% Print data019
Null_datafile = '/Volumes/Analysis/2015-02-24-2/dNull-norefit_nps/data019-from-data012_data014_data019_data006_data010_data016';
neuronPath = [Null_datafile,sprintf('/data019-from-data012_data014_data019_data006_data010_data016.neurons')];


ref_cell_number=1;
nConditions=3;
condDuration=12;

[spkColl,spkCondColl019,h]=plot_raster_script_pc2015_02_24_2_light(InterestingCell_vis_id(ref_cell_number),nConditions,condDuration,cond_str,neuronPath)

% Print data014
Null_datafile = '/Volumes/Analysis/2015-02-24-2/dNull-norefit_nps/data014-from-data012_data014_data019_data006_data010_data016';
neuronPath = [Null_datafile,sprintf('/data014-from-data012_data014_data019_data006_data010_data016.neurons')];

imov=1;
ref_cell_number=1;
nConditions=3;
condDuration=12;


[spkColl,spkCondColl014,h]=plot_raster_script_pc2015_02_24_2_light(InterestingCell_vis_id(ref_cell_number),nConditions,condDuration,cond_str,neuronPath)

% Print data012
Null_datafile = '/Volumes/Analysis/2015-02-24-2/dNull-norefit_nps/data012-from-data012_data014_data019_data006_data010_data016';
neuronPath = [Null_datafile,sprintf('/data012-from-data012_data014_data019_data006_data010_data016.neurons')];

imov=1;
ref_cell_number=1;
nConditions=3;
condDuration=12;

[spkColl,spkCondColl012,h]=plot_raster_script_pc2015_02_24_2_light(InterestingCell_vis_id(ref_cell_number),nConditions,condDuration,cond_str,neuronPath)

col='rkmrkmrkm';
databaseShift = 3*30;
nTrials1=30;
h=figure('Color','w');
cnt=0;

for icond=1:3
    xPoints=spkCondColl019(icond).xPoints;
    yPoints=spkCondColl019(icond).yPoints;
plot(xPoints*1/20000, yPoints+cnt,'r');
cnt=cnt+max(yPoints);
hold on

    xPoints=spkCondColl014(icond).xPoints;
    yPoints=spkCondColl014(icond).yPoints;
plot(xPoints*1/20000, yPoints+cnt,'k');
cnt=cnt+max(yPoints);
hold on

    xPoints=spkCondColl012(icond).xPoints;
    yPoints=spkCondColl012(icond).yPoints;
plot(xPoints*1/20000, yPoints+cnt,'m');
cnt=cnt+max(yPoints);
hold on
end

xlim([0,condDuration]);
ylim([0,cnt]);
title('data012,014,019 - order: Null2,Null1,Original');

   if(~isdir(sprintf('/Volumes/Analysis/nishal/analyse_2015_02_24_2/data006/NullCells1')))
    mkdir(sprintf('/Volumes/Analysis/nishal/analyse_2015_02_24_2/data006/NullCells1'));
    end
   s=hgexport('readstyle','ras5');
   hgexport(h,sprintf('/Volumes/Analysis/nishal/analyse_2015_02_24_2/data006/NullCells1/CellID_%d.eps',cellListOrig(icell)),s);
  

end

%% Joint analysis - from 010

WN_datafile = '2015-02-24-2/streamed/data010/data010';
datarun1=load_data(WN_datafile);
datarun1=load_params(datarun1);
datarun1=load_ei(datarun1,datarun1.cell_ids);

WN_datafile_concat = '2015-02-24-2/dNull-norefit_nps/data010-from-data012_data014_data019_data006_data010_data016/data010-from-data012_data014_data019_data006_data010_data016';
datarun2=load_data(WN_datafile_concat);
datarun2=load_params(datarun2);
datarun2.names.rrs_ei_path='/Volumes/Analysis/2015-02-24-2/dNull-norefit_nps/dNull-norefit_nps.ei';
datarun2=load_ei(datarun2,datarun2.cell_ids);

[cell_list_map, failed_cells] =map_ei(datarun1,datarun2);

NullCells1=[5272,3856,271,7067,5963,2528];
NullCells2 = [5011,1458,1832,7126,2596];

NullCells1_2=[];
NullCells1_1=[];
icnt=1;
for iList=1:length(NullCells1)
x=cell_list_map(datarun1.cell_ids==NullCells1(iList));
try
NullCells1_2(icnt)=x{1};
NullCells1_1(icnt)=NullCells1(icnt);
icnt=icnt+1;
catch
display('No match');
end

end

NullCells2_2=[];
NullCells2_1=[];
icnt=1;
for iList=1:length(NullCells2)
x=cell_list_map(datarun1.cell_ids==NullCells2(iList));
try
NullCells2_2(icnt)=x{1};
NullCells2_1(icnt)=NullCells2(icnt);
icnt=icnt+1;
catch
display('No match');    
end
end

%
cellListOrig=NullCells2_1;
cellList=NullCells2_2;
for icell=1:length(cellList);
    InterestingCell_vis_id=cellList(icell);


% Print data019
Null_datafile = '/Volumes/Analysis/2015-02-24-2/dNull-norefit_nps/data019-from-data012_data014_data019_data006_data010_data016';
neuronPath = [Null_datafile,sprintf('/data019-from-data012_data014_data019_data006_data010_data016.neurons')];


ref_cell_number=1;
nConditions=3;
condDuration=12;

[spkColl,spkCondColl019,h]=plot_raster_script_pc2015_02_24_2_light(InterestingCell_vis_id(ref_cell_number),nConditions,condDuration,cond_str,neuronPath)

% Print data014
Null_datafile = '/Volumes/Analysis/2015-02-24-2/dNull-norefit_nps/data014-from-data012_data014_data019_data006_data010_data016';
neuronPath = [Null_datafile,sprintf('/data014-from-data012_data014_data019_data006_data010_data016.neurons')];

imov=1;
ref_cell_number=1;
nConditions=3;
condDuration=12;


[spkColl,spkCondColl014,h]=plot_raster_script_pc2015_02_24_2_light(InterestingCell_vis_id(ref_cell_number),nConditions,condDuration,cond_str,neuronPath)

% Print data012
Null_datafile = '/Volumes/Analysis/2015-02-24-2/dNull-norefit_nps/data012-from-data012_data014_data019_data006_data010_data016';
neuronPath = [Null_datafile,sprintf('/data012-from-data012_data014_data019_data006_data010_data016.neurons')];

imov=1;
ref_cell_number=1;
nConditions=3;
condDuration=12;

[spkColl,spkCondColl012,h]=plot_raster_script_pc2015_02_24_2_light(InterestingCell_vis_id(ref_cell_number),nConditions,condDuration,cond_str,neuronPath)

col='rkmrkmrkm';
databaseShift = 3*30;
nTrials1=30;
h=figure('Color','w');
cnt=0;

for icond=1:3
    xPoints=spkCondColl019(icond).xPoints;
    yPoints=spkCondColl019(icond).yPoints;
plot(xPoints*1/20000, yPoints+cnt,'r');
cnt=cnt+max(yPoints);
hold on

    xPoints=spkCondColl014(icond).xPoints;
    yPoints=spkCondColl014(icond).yPoints;
plot(xPoints*1/20000, yPoints+cnt,'k');
cnt=cnt+max(yPoints);
hold on

    xPoints=spkCondColl012(icond).xPoints;
    yPoints=spkCondColl012(icond).yPoints;
plot(xPoints*1/20000, yPoints+cnt,'m');
cnt=cnt+max(yPoints);
hold on
end

xlim([0,condDuration]);
ylim([0,cnt]);
title('data012,014,019 - order: Null2,Null1,Original');

   if(~isdir(sprintf('/Volumes/Analysis/nishal/analyse_2015_02_24_2/data010/NullCells2')))
    mkdir(sprintf('/Volumes/Analysis/nishal/analyse_2015_02_24_2/data010/NullCells2'));
    end
   s=hgexport('readstyle','ras5');
   hgexport(h,sprintf('/Volumes/Analysis/nishal/analyse_2015_02_24_2/data010/NullCells2/CellID_%d.eps',cellListOrig(icell)),s);
  

end