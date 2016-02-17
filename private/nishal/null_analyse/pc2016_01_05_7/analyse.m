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

interestingConditions=[1,2,3,4,5];

dataRuns_ON_additivity = [3,4,5,6,7];
movies_ON_additivity = [1,2,3,4,6];


%% additivity experiment -  ON

dataRuns =dataRuns_ON_additivity;

WN_datafile = '/Volumes/Analysis/2016-01-05-7/d00_20-norefit/data000/data000';


datarun=load_data(WN_datafile)
datarun=load_params(datarun)

cellTypeId = 1
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
    Null_datafile = sprintf('/Volumes/Analysis/2016-01-05-7/d00_20-norefit/data0%02d',dataRuns(idata));
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
    
    if(~isdir(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2016_01_05_7/d_add_00/CellType_%s',datarun.cell_types{cellTypeId}.name)))
        mkdir(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2016_01_05_7/d_add_00/CellType_%s',datarun.cell_types{cellTypeId}.name));
    end
   s=hgexport('readstyle','raster');
   hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2016_01_05_7/d_add_00/CellType_%s/CellID_%d.eps',datarun.cell_types{cellTypeId}.name,InterestingCell_vis_id(ref_cell_number)),s);
end


%% Analyze runs from data003 - Contrast controlled Nulling
% Condition strings
nConditions=8;
condDuration=1200/4;
cond_str=cell(nConditions,1);
cond_str{1}='Original';
cond_str{2}='Null ';

interestingConditions=[1,2,3,4,5,6,7,8,9,10,11,12];



% contrast controlled nulling [0.36,0.24,0.12] WN and WN nulled
dataRuns_ON_additivity = [8,9,10,14,15,16,11,12,13,17,18,19];
movies_ON_additivity = [3,4,6,15,16,18,11,12,14,19,20,22];
% dataRuns_ON_additivity = [8,9,10,14,15,16];
% movies_ON_additivity = [3,4,6,15,16,18];
selected_cells = [752,1161,2063,2313,3320,4145,4337,4941];
%% additivity experiment -  ON

dataRuns =dataRuns_ON_additivity;

WN_datafile = '/Volumes/Analysis/2016-01-05-7/d00_20-norefit/data000/data000';

datarun=load_data(WN_datafile)
datarun=load_params(datarun)

cellTypeId = 1
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
    Null_datafile = sprintf('/Volumes/Analysis/2016-01-05-7/d00_20-norefit/data0%02d',dataRuns(idata));
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
    
    if(~isdir(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2016_01_05_7/d_contrast_null_00_b/CellType_%s',datarun.cell_types{cellTypeId}.name)))
        mkdir(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2016_01_05_7/d_contrast_null_00_b/CellType_%s',datarun.cell_types{cellTypeId}.name));
    end
   s=hgexport('readstyle','raster');
   hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2016_01_05_7/d_contrast_null_00_b/CellType_%s/CellID_%d.eps',datarun.cell_types{cellTypeId}.name,InterestingCell_vis_id(ref_cell_number)),s);
end
