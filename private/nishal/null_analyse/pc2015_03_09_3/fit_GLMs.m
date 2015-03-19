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


WN_datafile = '2015-03-09-3/streamed/data008/data008';
WN_datafile_short='2015-03-09-3/streamed/data008/data008';
movie_xml = 'RGB-8-4-0.48-11111';
stim_len=900;% in seconds

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
    
    
       if(~isdir(sprintf('/Volumes/Analysis/nishal/analyse_2015_03_09_3/data008/CellType_%s',datarun.cell_types{cellTypeId}.name)))
    mkdir(sprintf('/Volumes/Analysis/nishal/analyse_2015_03_09_3/data008/CellType_%s',datarun.cell_types{cellTypeId}.name));
    end
 
  
  xx=datarun.cell_types{cellTypeId}.name;
  xx(xx==' ') = '';
  dsave=sprintf('/Volumes/Analysis/nishal/analyse_2015_03_09_3/data008/CellType_%s',xx);

 
   fittedGLM=glm_fit_from_WNrun({cellID}, WN_datafile_short, movie_xml, stim_len,dsave);
   InterestingCell_vis_id(ref_cell_number)

end


end

%% Get a cell and fit GLM


WN_datafile = '2015-03-09-3/streamed/data000/data000';
WN_datafile_short='2015-03-09-3/streamed/data000/data000';
movie_xml = 'RGB-8-1-0.48-11111';
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
    
      if(~isdir(sprintf('/Volumes/Analysis/nishal/analyse_2015_03_09_3/data000/CellType_%s',datarun.cell_types{cellTypeId}.name)))
    mkdir(sprintf('/Volumes/Analysis/nishal/analyse_2015_03_09_3/data000/CellType_%s',datarun.cell_types{cellTypeId}.name));
    end
 
  xx=datarun.cell_types{cellTypeId}.name;
  xx(xx==' ') = '';
  dsave=sprintf('/Volumes/Analysis/nishal/analyse_2015_03_09_3/data000/CellType_%s',xx);
 
   fittedGLM=glm_fit_from_WNrun({cellID}, WN_datafile_short, movie_xml, stim_len,dsave);
   InterestingCell_vis_id(ref_cell_number)
  
end


end


