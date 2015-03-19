%% fit GLMs

%% add path to nora's folder for GLM code
location_of_git_repo='/home/vision/Nishal/matlab/';
addpath(genpath([location_of_git_repo '/private/nora']));

% Java library
javaaddpath('/Volumes/Lab/Development/vision7/Vision.app/Contents/Resources/Java/Vision.jar');
addpath(genpath('~/Nishal/matlab/private/nishal/create_act_2/'));
addpath(genpath('~/Nishal/matlab/private/nishal/Test suite/GLM/'));
addpath(genpath('~/Nishal/matlab/private/nishal/Test suite/GLM2/'));
addpath(genpath('~/Nishal/matlab/code'));

%% If want to do run code for experiment day
% rmpath(genpath('~/Nishal/matlab/private/nishal/create_act_2/'));
% addpath(genpath('~/Nishal/matlab/private/nishal/create_act_2_exp_repo/create_act_2_pc2014_11_05_2/'));
% code_version='old';

%% Get a cell and fit GLM

% 
% WN_datafile = '2015-03-09-2/streamed/data038/data038';
% WN_datafile_short='2015-03-09-2/streamed/data038/data038';
% movie_xml = 'BW-8-2-0.48-11111-40x40';
% stim_len=1800;% in seconds
% 
% datarun=load_data(WN_datafile)
% datarun=load_params(datarun)
% 
% 
% for cellTypeId=1:2; % 1 for On Parasols, 2 for Off parasols
% InterestingCell_vis_id=[];
% for icellType=cellTypeId
%     icellType 
%     InterestingCell_vis_id=[InterestingCell_vis_id,datarun.cell_types{icellType}.cell_ids];
% end 
% 
% 
% 
% 
% for ref_cell_number=1:length(InterestingCell_vis_id); %11
%     close all
%     cellID=InterestingCell_vis_id(ref_cell_number);
%    fittedGLM=glm_fit_from_WNrun({cellID}, WN_datafile_short, movie_xml, stim_len);
%    InterestingCell_vis_id(ref_cell_number)
%     if(~isdir(sprintf('/Volumes/Analysis/nishal/analyse_2015_03_09_2/data038/CellType_%s',datarun.cell_types{cellTypeId}.name)))
%     mkdir(sprintf('/Volumes/Analysis/nishal/analyse_2015_03_09_2/data038/CellType_%s',datarun.cell_types{cellTypeId}.name));
%     end
%  
%    save(sprintf('/Volumes/Analysis/nishal/analyse_2015_03_09_2/data038/CellType_%s/CellID_%d.mat',datarun.cell_types{cellTypeId}.name,InterestingCell_vis_id(ref_cell_number)));
%  
% end
% 
% 
% end


%% Get a cell and fit GLM


WN_datafile = '2015-03-09-2/data031/data031';
WN_datafile_short='2015-03-09-2/data031/data031';
movie_xml = 'BW-8-6-0.48-11111-40x40';
stim_len=1800;% in seconds

datarun=load_data(WN_datafile)
datarun=load_params(datarun)


for cellTypeId=1:2; % 1 for On Parasols, 2 for Off parasols
InterestingCell_vis_id=[];
for icellType=cellTypeId
    icellType 
    InterestingCell_vis_id=[InterestingCell_vis_id,datarun.cell_types{icellType}.cell_ids];
end 




for ref_cell_number=1:length(InterestingCell_vis_id); %11
    close all
    cellID=InterestingCell_vis_id(ref_cell_number);
   fittedGLM=glm_fit_from_WNrun({cellID}, WN_datafile_short, movie_xml, stim_len);
   InterestingCell_vis_id(ref_cell_number)
    if(~isdir(sprintf('/Volumes/Analysis/nishal/analyse_2015_03_09_2/data031/CellType_%s',datarun.cell_types{cellTypeId}.name)))
    mkdir(sprintf('/Volumes/Analysis/nishal/analyse_2015_03_09_2/data031/CellType_%s',datarun.cell_types{cellTypeId}.name));
    end
 
   save(sprintf('/Volumes/Analysis/nishal/analyse_2015_03_09_2/data031/CellType_%s/CellID_%d.mat',datarun.cell_types{cellTypeId}.name,InterestingCell_vis_id(ref_cell_number)));
 
end


end