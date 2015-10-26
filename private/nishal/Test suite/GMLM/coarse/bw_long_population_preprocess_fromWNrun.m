
location_of_git_repo='/home/vision/Nishal/matlab/';
addpath(genpath([location_of_git_repo '/private/nora']));

% Java library

addpath(('/Volumes/Lab/Users/bhaishahster/GITs/matlab/private/nishal/create_act2'));
javaaddpath('/Volumes/Lab/Development/vision7/Vision.app/Contents/Resources/Java/Vision.jar');
addpath(genpath('/Volumes/Lab/Users/bhaishahster/GITs/matlab/private/nishal/Test suite/GLM/'));
addpath(genpath('/Volumes/Lab/Users/bhaishahster/GITs/matlab/private/nishal/Test suite/GLM2/'));
addpath(genpath('/Volumes/Lab/Users/bhaishahster/GITs/matlab/private/nishal'));
addpath(genpath('/Volumes/Lab/Users/bhaishahster/GITs/matlab/code'));
addpath(genpath('/Volumes/Lab/Users/bhaishahster/GITs/matlab/private/nishal/create_act_2/'));
%% Dataset details

WN_datafile = '2005-04-26-2/data002/data002';
WN_datafile_short=WN_datafile;
movie_xml = 'RGB-20-1-0.48-22222';
stim_length=1800;% 
stix=20;
%% Load stimuli


datarun=load_data(WN_datafile)
datarun=load_params(datarun)
datarun=load_sta(datarun);
datarun=load_neurons(datarun);

icell_l=0;
cells = [datarun.cell_types{1}.cell_ids];
ttf_log = zeros(30,length(cells));
totalMaskAccept_log = zeros((640/stix)*(320/stix),length(cells));
Y = zeros(length(cells),216000);

for cellID = cells
    icell_l=icell_l+1
extract_movie_response4;
ttf_log(:,icell_l) =ttf;
totalMaskAccept_log(:,icell_l) = totalMaskAccept(:);
Y(icell_l,:) = spksGen;
end

 ttf_avg = mean(ttf_log,2);

 mov=squeeze(mean(mov,3));
 maskedMovdd= filterMov(mov,ones(size(mov,1),size(mov,2)),squeeze(ttf_avg));

  save('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2005_04_26_2/ON_parasol.mat','cells','maskedMovdd','Y','ttf_log','ttf_avg','totalMaskAccept_log','mov','-v7.3')