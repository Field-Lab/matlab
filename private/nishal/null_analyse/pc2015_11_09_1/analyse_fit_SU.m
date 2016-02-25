
%% 
WN_datafile = '/Volumes/Analysis/2015-11-09-1/d14_36-norefit/data009/data009';
datarun=load_data(WN_datafile)
datarun=load_params(datarun)


cellTypeId = 1%OFF:2 .. ON:1.. 
cellID_select= datarun.cell_types{cellTypeId}.cell_ids; % 51 ? [6741,2176,6961]


 %% fit ASM - from data002 (low contrast)
 
 % OFF
cellTypeId = 2%OFF:2 .. ON:1.. 
cellID_select= datarun.cell_types{cellTypeId}.cell_ids; % 51 ? [6741,2176,6961]

WN_datafile = '/Volumes/Analysis/2015-11-09-1/d14_36-norefit/data026/data026'
movie_xml = 'BW-8-4-0.48-11111';
stim_length=900;% 
cellID_list = [cellID_select];%1006,1411,2086,2341,2686,4096,4276,4831,5371,5566,5596,5866,7156,7216,7381,7490];
nSU_list= [1,2,3,4,5,7,10];
destination= 'pc2015_11_09_1_analysis_fits/SUs_data026/OFF '
sta_depth=30; % CAREFUL TO CHANGE IT TO LONG .. FOR SLOWER CELLS
contrast_factor=0.5;
[~] = get_gmlm_sta2(WN_datafile,movie_xml,stim_length,cellID_list,nSU_list,destination,sta_depth,contrast_factor)
save(['/Volumes/Lab/Users/bhaishahster/',destination,'/sus.mat'],'stas_t','stas_r');



 
 % ON
cellTypeId = 1%OFF:2 .. ON:1.. 
cellID_select= datarun.cell_types{cellTypeId}.cell_ids; % 51 ? [6741,2176,6961]

WN_datafile = '/Volumes/Analysis/2015-11-09-1/d14_36-norefit/data026/data026'
movie_xml = 'BW-8-4-0.48-11111';
stim_length=900;% 
cellID_list = [cellID_select];%1006,1411,2086,2341,2686,4096,4276,4831,5371,5566,5596,5866,7156,7216,7381,7490];
nSU_list= [1,2,3,4,5,7,10];
destination= 'pc2015_11_09_1_analysis_fits/SUs_data026/ON'
sta_depth=30; % CAREFUL TO CHANGE IT TO LONG .. FOR SLOWER CELLS
contrast_factor=0.5;
[~] = get_gmlm_sta2(WN_datafile,movie_xml,stim_length,cellID_list,nSU_list,destination,sta_depth,contrast_factor)