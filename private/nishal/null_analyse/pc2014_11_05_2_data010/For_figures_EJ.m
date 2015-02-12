
WN_datafile = '2014-11-05-2/data009/data009';
Null_datafile = '/Volumes/Analysis/2014-11-05-2/data010';
%WN_mov='/Volumes/Analysis/stimuli/white-noise-xml/RGB-10-2-0.48-11111.xml';
WN_datafile_full = '/Volumes/Analysis/2014-11-05-2/data009';
imov=10;


datarun=load_data(WN_datafile)
datarun=load_params(datarun)


% cellTypeId=[2]; % 1 for On Parasols, 2 for Off parasols
% InterestingCell_vis_id=[];
% for icellType=cellTypeId
%     icellType 
%     InterestingCell_vis_id=[InterestingCell_vis_id,datarun.cell_types{icellType}.cell_ids];
% end 
% cellTypeUsed=cellTypeId*ones(length(InterestingCell_vis_id),1);
InterestingCell_vis_id = [275,1382,2181,2255,2389,2716,2747,2957,3362,4567,5056,5057,5928,5929,6543,7276,7742]
for ref_cell_number=[16]%1:length(InterestingCell_vis_id); %11
   
    close all
    PltCnd=[2,1];
    [spkColl,spkCondColl]=plot_raster_script_withGLM_pred(datarun,WN_datafile,WN_datafile_full,Null_datafile,InterestingCell_vis_id,imov,ref_cell_number,nConditions,condDuration,cond_str,PltCnd);
    ref_cell_number
    CondtoWrite=[1,2];
    makeConditionMovies( cond_str,condMovies,datarun,InterestingCell_vis_id(ref_cell_number),CondtoWrite);
    
   % plot_mosaic(datarun,InterestingCell_vis_id,ref_cell_number)
% pause
end