
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

WN_datafile = '2015-10-29-2/d00_36-norefit/data002/data002'; % data002 RGB-8-2-0.24-11111, 30 min, all streamed , low contrast (0.24)
WN_datafile_short=WN_datafile;
movie_xml = 'RGB-8-2-0.48-11111';
stim_length=1800;% 
%% Load stimuli


datarun=load_data(WN_datafile)
datarun=load_params(datarun)
datarun=load_sta(datarun);
datarun=load_neurons(datarun);

icell_l=0;
cells = [datarun.cell_types{2}.cell_ids];
ttf_log = zeros(30,length(cells));
totalMaskAccept_log = zeros(80*40,length(cells));
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

  save('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/data002_analysis/Off_parasol.mat','maskedMovdd','Y','ttf_log','ttf_avg','totalMaskAccept_log','-v7.3')
  
  
  %% fit ASM
  
path = '/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/data002_analysis'
load([path,'/Off_parasol.mat']);
binnedSpikeResponses_coll = Y;
total_mask_log = totalMaskAccept_log;

cells = double(cells);
cellsChoose = (cells ==3287) | (cells ==3318 ) | (cells ==3155) | (cells ==3066);

cellsChoose=logical(cellsChoose);
mask = sum(total_mask_log(:,cellsChoose),2)~=0;
ttf = ttf_avg;


ifit=0;
for fitNnum=1:1
for Ns=[4,8,6,10,12,14,16];
%fitASM_pop = fitASM_EM_Population(1000*maskedMovdd,binnedSpikeResponses_coll(cellsChoose,:),Ns,mask);
ifit=ifit+1
lam_start = 0.1;
gamma1=0.0000;
gamma2 = 0;
initVal=[];
%initVal.K = maskedMovdd(mask,:)*binnedSpikeResponses_coll(cellsChoose,:)';
%initVal.B = diag(ones(Ns,1));
[fitASM_pop,fval] = fitASM_EM_Population_sparse_split_admm(maskedMovdd,binnedSpikeResponses_coll(cellsChoose,:),Ns,mask,gamma1,gamma2,lam_start,initVal);
B_use= plotSU_withcells(fitASM_pop.K,mask,total_mask_log(:,cellsChoose),exp(fitASM_pop.B));
save(sprintf('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/data002_sharedSU/Off_parasol_%d.mat',ifit),'Ns','fitASM_pop','fval','mask','gamma1','gamma2','lam_start','initVal','cellsChoose');
pause(0.5);
end
end 
  