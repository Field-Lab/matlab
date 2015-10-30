%% do streaming analysis of a new dataset

% Use Java on server (where data resides)
% Goto /Volumes/Lab/Development/vision7/
% java -jar Vision.jar
%% 

dataRuns = [3,4,17,21,19,23,25,29,27,31];

WN_datafile_full = '/Volumes/Analysis/2015-10-29-2/streamed/data001/data001';

datarun=load_data(WN_datafile_full)
datarun=load_params(datarun)

cellTypeId = 2%8,1;
InterestingCell_vis_id=7126%751%datarun.cell_types{cellTypeId}.cell_ids; %[556,1278,1384,407,1516,2150,2401,3361,4066,4683,5611,6106,6005,7246,7562,3946];
cellTypeUsed=cellTypeId*ones(length(InterestingCell_vis_id),1);

condDuration=10;
nConditions=1;

cols='rkrkrkrkrkrkkrkrkrkr';
spkCondColl=cell(8,1);

for ref_cell_number=1%1:length(InterestingCell_vis_id); %11
    close all
    cellID=InterestingCell_vis_id(ref_cell_number)
        
    for idata=1:length(dataRuns)
    Null_datafile = sprintf(sprintf('/Volumes/Analysis/2015-10-29-2/streamed/data0%02d-from-data001/',dataRuns(idata))); % depends on idata .. 
    neuronPath = [Null_datafile,sprintf('data0%02d-from-data001.neurons',dataRuns(idata))];
    [spkColl,spkCondColl{idata},h]=plot_raster_script_pc2015_09_23_0_light(cellID,nConditions,condDuration,cond_str,neuronPath);
    end
    
    
    h=figure;
    for irun = 1:length(dataRuns)
    plot(spkCondColl{irun}.xPoints/20000,spkCondColl{irun}.yPoints - (irun-1)*30,cols(irun));
    hold on;
    end
    set(gca,'yTick',[]);
    ylim([-(length(dataRuns)-1)*30,30]);
    InterestingCell_vis_id(ref_cell_number)
    
 %   if(~isdir(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_10_06_0/dsu_null_29_31_33_35/CellType_%s',datarun.cell_types{cellTypeId}.name)))
 %      mkdir(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_10_06_0/dsu_null_29_31_33_35/CellType_%s',datarun.cell_types{cellTypeId}.name));
 %   end
 %  s=hgexport('readstyle','raster');
 %  hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_10_06_0/dsu_null_29_31_33_35/CellType_%s/CellID_%d.eps',datarun.cell_types{cellTypeId}.name,InterestingCell_vis_id(ref_cell_number)),s);
end





