 addpath(genpath('../null_analyse/'));
addpath(genpath('../null_analyse/analyse_functions'));
%startup_null_analyse_tenessee
startup_null_analyse_bertha

%% data003 from data 000
% Condition strings

cond_str=cell(6,1);
cond_str{1}='Original';
cond_str{2}='Null ';

%

WN_datafile = '2015-05-27-10/streamed/data000/data000';
Null_datafile = '/Volumes/Analysis/2015-05-27-10/data006_from_data000_s_nps';
%WN_mov='/Volumes/Analysis/stimuli/white-noise-xml/RGB-10-2-0.48-11111.xml';
WN_datafile_full = '2015-05-27-10/streamed/data000/data000';
neuronPath = [Null_datafile,sprintf('/data006_from_data000_s_nps.neurons')];
imov=2;

datarun=load_data(WN_datafile)
datarun=load_params(datarun)


cellTypeId=[1]; % 1 for On Parasols, 2 for Off parasols
InterestingCell_vis_id=[];
for icellType=cellTypeId
    icellType
    InterestingCell_vis_id=[InterestingCell_vis_id,datarun.cell_types{icellType}.cell_ids];
end

cellTypeUsed=cellTypeId*ones(length(InterestingCell_vis_id),1);
NullCells1=InterestingCell_vis_id; 
NullCells2=[];

condDuration=10;
nConditions=3;
for ref_cell_number=1:length(InterestingCell_vis_id); %11
    close all
    cellID=InterestingCell_vis_id(ref_cell_number);
     [spkColl,spkCondColl,h]=plot_raster_script_pc2015_05_27_10_light(cellID,nConditions,condDuration,cond_str,neuronPath);
    set(gca,'YTick',[]);
    plot_mosaic_pc2015_03_09_2(datarun,InterestingCell_vis_id,ref_cell_number,NullCells1,NullCells2);
   
    InterestingCell_vis_id(ref_cell_number)
%     
%     if(~isdir(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_05_27_10/data003/CellType_%s/CellID_%d',datarun.cell_types{cellTypeId}.name,InterestingCell_vis_id(ref_cell_number))))
%         mkdir(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_05_27_10/data003/CellType_%s/CellID_%d',datarun.cell_types{cellTypeId}.name,InterestingCell_vis_id(ref_cell_number)));
%     end
    if(~isdir(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_05_27_10/data006/CellType_%s',datarun.cell_types{cellTypeId}.name)))
        mkdir(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_05_27_10/data006/CellType_%s',datarun.cell_types{cellTypeId}.name));
    end
   s=hgexport('readstyle','ras_mos4');
   hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_05_27_10/data006/CellType_%s/CellID_%d.eps',datarun.cell_types{cellTypeId}.name,InterestingCell_vis_id(ref_cell_number)),s);
   % pause
end

