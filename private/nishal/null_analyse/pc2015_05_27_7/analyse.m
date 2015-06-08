 addpath(genpath('../null_analyse/'));
addpath(genpath('../null_analyse/analyse_functions'));
%startup_null_analyse_tenessee
startup_null_analyse_bertha

%% data004 from data 002
% Condition strings

cond_str=cell(6,1);
cond_str{1}='Original';
cond_str{2}='Null ';
cond_str{3}='';
cond_str{4}='';
cond_str{5}='';
cond_str{6}='';    
%

WN_datafile = '2015-05-27-7/streamed/data002/data002';
Null_datafile = '/Volumes/Analysis/2015-05-27-7/data004_from_data002_s_nps';
%WN_mov='/Volumes/Analysis/stimuli/white-noise-xml/RGB-10-2-0.48-11111.xml';
WN_datafile_full = '2015-05-27-7/streamed/data002/data002';
neuronPath = [Null_datafile,sprintf('/data004_from_data002_s_nps.neurons')];
imov=2;

datarun=load_data(WN_datafile)
datarun=load_params(datarun)


cellTypeId=[1]; % 1 for On Parasols, 2 for Off parasols
InterestingCell_vis_id=[];
for icellType=cellTypeId
    icellType
    InterestingCell_vis_id=[InterestingCell_vis_id,datarun.cell_types{icellType}.cell_ids];
end
InterestingCell_vis_id = [1561,1921,2821,3826,4141,4816,4922,5761,5806,6091,7563,7426,751,1486,1681,3229,4081,4381,5072,6661];
cellTypeUsed=cellTypeId*ones(length(InterestingCell_vis_id),1);
NullCells1=InterestingCell_vis_id; 
NullCells2=[];

condDuration=10;
nConditions=6;
for ref_cell_number=1:length(InterestingCell_vis_id); %11
    close all
    cellID=InterestingCell_vis_id(ref_cell_number);
     [spkColl,spkCondColl,h]=plot_raster_script_pc2015_05_27_7_light(cellID,nConditions,condDuration,cond_str,neuronPath);
    set(gca,'YTick',[]);
   % plot_mosaic_pc2015_05_27_7(datarun,InterestingCell_vis_id,ref_cell_number,NullCells1,NullCells2);
   
    InterestingCell_vis_id(ref_cell_number)
%     
%     if(~isdir(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_05_27_10/data003/CellType_%s/CellID_%d',datarun.cell_types{cellTypeId}.name,InterestingCell_vis_id(ref_cell_number))))
%         mkdir(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_05_27_10/data003/CellType_%s/CellID_%d',datarun.cell_types{cellTypeId}.name,InterestingCell_vis_id(ref_cell_number)));
%     end
    if(~isdir(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_05_27_7/data004/CellType_%s',datarun.cell_types{cellTypeId}.name)))
        mkdir(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_05_27_7/data004/CellType_%s',datarun.cell_types{cellTypeId}.name));
    end
   s=hgexport('readstyle','ras_mos4');
   hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_05_27_7/data004/CellType_%s/CellID_%d.eps',datarun.cell_types{cellTypeId}.name,InterestingCell_vis_id(ref_cell_number)),s);
    pause
end

