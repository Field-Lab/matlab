
WN_datafile = '2015-03-09-2/streamed/data038/data038';
Null_datafile = '/Volumes/Analysis/2015-03-09-2/data042-from-data038_streamed_nps';
%WN_mov='/Volumes/Analysis/stimuli/white-noise-xml/RGB-10-2-0.48-11111.xml';
WN_datafile_full = '/Volumes/Analysis/2015-03-09-2/data038/data038';
neuronPath = [Null_datafile,sprintf('/data042-from-data038_streamed_nps.neurons')];
imov=10;



datarun=load_data(WN_datafile)
datarun=load_params(datarun);

InterestingCell_vis_id = [1531]
for ref_cell_number=[1]
   
    close all
    PltCnd=[5,1];
    [spkColl,spkCondColl]=plot_raster_script_withGLM_pred(datarun,WN_datafile,WN_datafile_full,Null_datafile,InterestingCell_vis_id,imov,ref_cell_number,nConditions,condDuration,cond_str,PltCnd);
    ref_cell_number
    CondtoWrite=[1,2];
    makeConditionMovies_pc2015_03_09_2( cond_str,condMov,datarun,InterestingCell_vis_id(ref_cell_number),CondtoWrite);
    
   % plot_mosaic(datarun,InterestingCell_vis_id,ref_cell_number)
% pause
end