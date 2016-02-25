addpath(genpath('../null_analyse/'));
addpath(genpath('../null_analyse/analyse_functions'));
%startup_null_analyse_tenessee
startup_null_analyse_bertha

%% 
WN_datafile = '/Volumes/Analysis/2016-02-17-1/d28_51-norefit/data028/data028';
datarun=load_data(WN_datafile)
datarun=load_params(datarun)



 %% fit ASM - from data002 (low contrast)
 
 % OFF
cellTypeId = 2%OFF:2 .. ON:1.. 
cellID_select= datarun.cell_types{cellTypeId}.cell_ids; % 51 ? [6741,2176,6961]

WN_datafile = '/Volumes/Analysis/2016-02-17-1/d28_51-norefit/data029/data029'
movie_xml = 'BW-8-4-0.48-11111-119.5';
stim_length=1700;% 
cellID_list = [cellID_select];%1006,1411,2086,2341,2686,4096,4276,4831,5371,5566,5596,5866,7156,7216,7381,7490];
nSU_list= [1,2,3,4,5,7,10];
destination= 'pc2016_02_17_1_analysis_fits/SUs_data029/OFF '
sta_depth=30; % CAREFUL TO CHANGE IT TO LONG .. FOR SLOWER CELLS
contrast_factor=0.5;
[~] = get_gmlm_sta2(WN_datafile,movie_xml,stim_length,cellID_list,nSU_list,destination,sta_depth,contrast_factor)
%save(['/Volumes/Lab/Users/bhaishahster/',destination,'/sus.mat'],'stas_t','stas_r');


!
 
 % ON
cellTypeId = 1%OFF:2 .. ON:1.. 
cellID_select= datarun.cell_types{cellTypeId}.cell_ids; % 51 ? [6741,2176,6961]

WN_datafile = '/Volumes/Analysis/2016-02-17-1/d28_51-norefit/data029/data029'
movie_xml = 'BW-8-4-0.48-11111-119.5';
stim_length=1700;% 
cellID_list = [cellID_select];%1006,1411,2086,2341,2686,4096,4276,4831,5371,5566,5596,5866,7156,7216,7381,7490];
cellID_list=cellID_list(cellID_list>5868);
nSU_list= [1,2,3,4,5,7,10];
destination= 'pc2016_02_17_1_analysis_fits/SUs_data029/ON'
sta_depth=30; % CAREFUL TO CHANGE IT TO LONG .. FOR SLOWER CELLS
contrast_factor=0.5;
[~] = get_gmlm_sta2(WN_datafile,movie_xml,stim_length,cellID_list,nSU_list,destination,sta_depth,contrast_factor)