addpath(genpath('../null_analyse/'));
addpath(genpath('../null_analyse/analyse_functions'));
%startup_null_analyse_tenessee
startup_null_analyse_bertha

%% Analyze runs from data001
% Condition strings
nConditions=8;
condDuration=1200/4;
cond_str=cell(nConditions,1);
cond_str{1}='Original';
cond_str{2}='Null ';

interestingConditions=[1,2,3,4];

dataRuns_OFF_additivity = [4,5,7,8];
dataRuns_ON_additivity = [4,6,7,9];
movies_OFF_additivity =[1,2,5,6];
movies_ON_additivity = [1,4,5,8];


%% additivity experiment -  OFF

dataRuns =dataRuns_OFF_additivity;

WN_datafile = '/Volumes/Analysis/2016-01-05-3/d00_37-norefit/data001/data001';


datarun=load_data(WN_datafile)
datarun=load_params(datarun)

cellTypeId = 2
InterestingCell_vis_id=datarun.cell_types{cellTypeId}.cell_ids; 
cellTypeUsed=cellTypeId*ones(length(InterestingCell_vis_id),1);

condDuration=10;
nConditions=1;

cols='rkrkrkrkrkrkkrkrkrkr';
spkCondColl=cell(length(dataRuns),1);

for ref_cell_number=1:length(InterestingCell_vis_id); %11
    close all
    cellID=InterestingCell_vis_id(ref_cell_number)
        
    for idata=1:length(dataRuns)
    Null_datafile = sprintf('/Volumes/Analysis/2016-01-05-3/d00_37-norefit/data0%02d',dataRuns(idata));
    neuronPath = [Null_datafile,sprintf('/data0%02d.neurons',dataRuns(idata))];
    [spkColl,spkCondColl{idata},h]=plot_raster_script_pc2016_01_05_3_light_photons(cellID,nConditions,condDuration,cond_str,neuronPath);
    end
    
    
    h=figure;
    for irun = 1:length(dataRuns)
    plot(spkCondColl{irun}.xPoints/20000,spkCondColl{irun}.yPoints - (irun-1)*30,cols(irun));
    hold on;
    end
    set(gca,'yTick',[]);
    ylim([-(length(dataRuns)-1)*30,30]);
    InterestingCell_vis_id(ref_cell_number)
    
    if(~isdir(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2016_01_05_3/d_add_01/CellType_%s',datarun.cell_types{cellTypeId}.name)))
        mkdir(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2016_01_05_3/d_add_01/CellType_%s',datarun.cell_types{cellTypeId}.name));
    end
   s=hgexport('readstyle','raster');
   hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2016_01_05_3/d_add_01/CellType_%s/CellID_%d.eps',datarun.cell_types{cellTypeId}.name,InterestingCell_vis_id(ref_cell_number)),s);
end


%% Analyze runs from data003 - Normal Additivity + Nulling
% Condition strings
nConditions=8;
condDuration=1200/4;
cond_str=cell(nConditions,1);
cond_str{1}='Original';
cond_str{2}='Null ';

interestingConditions=[1,2,3,4];

% Nulling, additivity experiments
% high contrast WN, high contrast Null, low contrast WN, low contrast Null,
% low contrast Null selected cells, WN+Null variance controlled, WN+Null
% Variance NOT controlled
dataRuns_OFF_additivity = [10,11,13,14,23,32,34,37]; 
dataRuns_ON_additivity = [10,12,13,15,24,33,35,37];
movies_OFF_additivity =[1,2,5,6,24,28,32];
movies_ON_additivity = [1,4,5,8,26,30,34];


% 
% % contrast controlled nulling [0.36,0.24,0.12] WN and WN nulled
% dataRuns_OFF_additivity = [16,17,19,21,25,26,37]; 
% dataRuns_ON_additivity = [16,18,19,22,25,27,37];
% movies_OFF_additivity =[9,10,13,14,17,18,1];
% movies_ON_additivity = [9,12,13,16,17,20,1];
%% additivity experiment -  OFF

dataRuns =dataRuns_ON_additivity;

WN_datafile = '/Volumes/Analysis/2016-01-05-3/d00_37-norefit/data003/data003'; % OFF Par : 14, ON Par : 1


datarun=load_data(WN_datafile)
datarun=load_params(datarun)

cellTypeId = 1%14
InterestingCell_vis_id=datarun.cell_types{cellTypeId}.cell_ids; 
cellTypeUsed=cellTypeId*ones(length(InterestingCell_vis_id),1);

condDuration=10;
nConditions=1;

cols='rkrkrkrkrkrkkrkrkrkr';
spkCondColl=cell(length(dataRuns),1);

for ref_cell_number=1:length(InterestingCell_vis_id); %11
    close all
    cellID=InterestingCell_vis_id(ref_cell_number)
        
    for idata=1:length(dataRuns)
    Null_datafile = sprintf('/Volumes/Analysis/2016-01-05-3/d00_37-norefit/data0%02d',dataRuns(idata));
    neuronPath = [Null_datafile,sprintf('/data0%02d.neurons',dataRuns(idata))];
    [spkColl,spkCondColl{idata},h]=plot_raster_script_pc2016_01_05_3_light_photons(cellID,nConditions,condDuration,cond_str,neuronPath);
    end
    
    
    h=figure;
    for irun = 1:length(dataRuns)
    plot(spkCondColl{irun}.xPoints/20000,spkCondColl{irun}.yPoints - (irun-1)*30,cols(irun));
    hold on;
    end
    set(gca,'yTick',[]);
    ylim([-(length(dataRuns)-1)*30,30]);
    InterestingCell_vis_id(ref_cell_number)
    
    if(~isdir(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2016_01_05_3/d_add_03/CellType_%s',datarun.cell_types{cellTypeId}.name)))
        mkdir(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2016_01_05_3/d_add_03/CellType_%s',datarun.cell_types{cellTypeId}.name));
    end
   s=hgexport('readstyle','raster');
   hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2016_01_05_3/d_add_03/CellType_%s/CellID_%d.eps',datarun.cell_types{cellTypeId}.name,InterestingCell_vis_id(ref_cell_number)),s);
end

%% Analyze runs from data003 - Contrast Controlled Nulling

% Condition strings
nConditions=8;
condDuration=1200/4;
cond_str=cell(nConditions,1);
cond_str{1}='Original';
cond_str{2}='Null ';

interestingConditions=[1,2,3,4];


% contrast controlled nulling [0.36,0.24,0.12] WN and WN nulled
dataRuns_OFF_additivity = [16,17,19,21,25,26,37]; 
dataRuns_ON_additivity = [16,18,19,22,25,27,37];
movies_OFF_additivity =[9,10,13,14,17,18,1];
movies_ON_additivity = [9,12,13,16,17,20,1];


%% additivity experiment -  OFF

dataRuns =dataRuns_ON_additivity;

WN_datafile = '/Volumes/Analysis/2016-01-05-3/d00_37-norefit/data003/data003'; % OFF Par : 14, ON Par : 1


datarun=load_data(WN_datafile)
datarun=load_params(datarun)

cellTypeId = 1%14
InterestingCell_vis_id=datarun.cell_types{cellTypeId}.cell_ids; 
cellTypeUsed=cellTypeId*ones(length(InterestingCell_vis_id),1);

condDuration=10;
nConditions=1;

cols='rkrkrkrkrkrkkrkrkrkr';
spkCondColl=cell(length(dataRuns),1);

for ref_cell_number=1:length(InterestingCell_vis_id); %11
    close all
    cellID=InterestingCell_vis_id(ref_cell_number)
        
    for idata=1:length(dataRuns)
    Null_datafile = sprintf('/Volumes/Analysis/2016-01-05-3/d00_37-norefit/data0%02d',dataRuns(idata));
    neuronPath = [Null_datafile,sprintf('/data0%02d.neurons',dataRuns(idata))];
    [spkColl,spkCondColl{idata},h]=plot_raster_script_pc2016_01_05_3_light_photons(cellID,nConditions,condDuration,cond_str,neuronPath);
    end
    
    
    h=figure;
    for irun = 1:length(dataRuns)
    plot(spkCondColl{irun}.xPoints/20000,spkCondColl{irun}.yPoints - (irun-1)*30,cols(irun));
    hold on;
    end
    set(gca,'yTick',[]);
    ylim([-(length(dataRuns)-1)*30,30]);
    InterestingCell_vis_id(ref_cell_number)
    
    if(~isdir(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2016_01_05_3/d_contrast_null_03/CellType_%s',datarun.cell_types{cellTypeId}.name)))
        mkdir(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2016_01_05_3/d_contrast_null_03/CellType_%s',datarun.cell_types{cellTypeId}.name));
    end
   s=hgexport('readstyle','raster');
   hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2016_01_05_3/d_contrast_null_03/CellType_%s/CellID_%d.eps',datarun.cell_types{cellTypeId}.name,InterestingCell_vis_id(ref_cell_number)),s);
end



