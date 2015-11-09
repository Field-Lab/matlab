
%% GMLM with MEL and EM for normal cells at coarser resolution .. 

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
%% Load data
WN_datafile = '2015-03-09-2/streamed/data038/data038';
WN_datafile_short='2015-03-09-2/streamed/data038/data038';
movie_xml = 'BW-8-2-0.48-11111-40x40';
stim_length=1800;% 


datarun=load_data(WN_datafile)
datarun=load_params(datarun)
datarun=load_sta(datarun);
datarun=load_neurons(datarun);

%% Find stimulus and response
cellID=1531;

user_STA_depth=30;
extract_movie_response2;
%% Load fits
data1 = load(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data038/Off parasol/detailed/CellID_%d/Cell%d_full.mat',cellID,cellID));
data2 = load(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data038/Off parasol/detailed/CellID_%d/Cell%d_full_su4_5.mat',cellID,cellID));

fitGMLM_log=cell(5,1);
for isu=1:3
fitGMLM_log{isu} = data1.fitGMLM_full2_log{isu};
end
for isu=4:5
fitGMLM_log{isu} = data2.fitGMLM_full2_log{isu};
end
mask=data2.totalMaskAccept2;

nSU=5;

fitGMLM = fitGMLM_log{nSU};

%% 
 mask=logical(mask);
K=zeros(sum(mask(:)),nSU);
B = ones(nSU,1);
for isu=1:nSU
K(:,isu) = fitGMLM.Linear.filter{isu}(1:end-1);
B(isu) = exp(fitGMLM.Linear.filter{isu}(end));
end


figure;
[B_use,h]= plotSU_withcells(K,logical(mask),logical(mask(:)),B);
suptitle('All SU');

sus = input('Which sus?');
figure;
isu=sus(1);
ss = zeros(40,40);ss(logical(mask))=K(:,isu);ss=ss';
plot_SU_with_idx(ss,mask);
title(sprintf('SU: %d',isu));

figure
jsu = sus(2);
ss = zeros(40,40);ss(logical(mask))=K(:,jsu);ss=ss';
plot_SU_with_idx(ss,mask);

title(sprintf('SU: %d',jsu));

Y_C = spksGen';mov = maskedMovdd;

%% Select movie segments where probe movie is strong

uProbe = zeros(sum(mask(:)),1);
cell_choose_num=[1];jcell=1;total_mask_log = mask(:);

suP = input('Testing which pixels in SU?')
q = -1*total_mask_log(:,cell_choose_num(jcell))*(sum(total_mask_log(:,cell_choose_num(jcell)))-length(suP)); uProbe = uProbe + q(mask);
uProbe(suP)=length(suP);
uProbe = uProbe/norm(uProbe);
 zz=zeros(40,40);zz(mask) = uProbe;figure;imagesc(zz');colormap gray;axis image
 
 % use uProbe to get interesting sections
 mov = maskedMovdd;
 
 movProbe = uProbe'*mov;
 thr = prctile(movProbe,90);
 probeTms = movProbe>thr;
 
 movProbeCos = movProbe./sqrt(sum(mov.^2,1));
 thr = prctile(movProbeCos,70); % Choose HOW? TODO!
 probeCosTms = movProbeCos>thr;
 
 probe = probeCosTms;
 sta = maskedMovdd(:,probe)*Y_C(jcell,probe)';
 xsta=zeros(40,40);xsta(mask)=sta;sta=xsta;

 figure;
 ssta = reshape(sta,[40,40]);
 total_mask_log=mask(:);cellsChoose=1;
[B_use,h]= plotSU_withcells(ssta(mask),mask,total_mask_log(:,cellsChoose),ones(nSU,1));
axis image
title('STA of stimuli chosen for analysis')




%% Do tests
 u1 = input('Test SU')%[8,9,18,19,30];
 u2 = input('Other pixels')%[10,20,32,31];
 upossible=input('Possible other pixels in SU')

 thperc=75;pts=7; type='geq';
 hh=figure;
 subplot(1,2,1);
 usu1 = u1;usu2=[];jcell=1;
 testSU_interaction(usu1,usu2,mask(:),mov,Y_C,jcell,probe,thperc,pts,type,'b')
 hold on;

 u1perm = perms(u1);u2perm=perms(u2); nperms=size(u2perm,1);npix = length(u1);
 for iperm=1:nperms
 usu1 = u1;
 usu2 = u2perm(iperm,1:(npix));
 testSU_interaction(usu1,usu2,mask(:),mov,Y_C,jcell,probe,thperc,pts,type,'r')
 hold on;
 end
 
 usu1 = u1;
 usu2 = upossible;
 testSU_interaction(usu1,usu2,mask(:),mov,Y_C,jcell,probe,thperc,pts,type,'k')
 hold on;
 
 
 subplot(1,2,2);
 thperc=25;pts=7; type='leq';
 usu1 = u1;usu2=[];
 testSU_interaction(usu1,usu2,mask(:),mov,Y_C,jcell,probe,thperc,pts,type,'b')
 hold on;

 u1perm = perms(u1);u2perm=perms(u2); nperms=size(u2perm,1);npix = length(u1);
 for iperm=1:nperms
     
 usu1 = u1;
 usu2 = u2perm(iperm,1:(npix));
 testSU_interaction(usu1,usu2,mask(:),mov,Y_C,jcell,probe,thperc,pts,type,'r')
 hold on;
 end
 
 usu1 = u1;
 usu2 = upossible;
 testSU_interaction(usu1,usu2,mask(:),mov,Y_C,jcell,probe,thperc,pts,type,'k')
 hold on;



 
 %% compute output non-linearity to ASM.
 
 figure;
 for nSU=1:5
     fitGMLM = fitGMLM_log{nSU};
 nTrials=1;
 [resp,lam] = predictGMLM_full(fitGMLM,maskedMovdd,nTrials);
 resp=resp';lam=lam';lam=lam/1200;
%  
% figure;
% scatter(lam,spksGen_hr,1*ones(length(spksGen_hr),1),'filled');hold on;

th_list=[];
for iprc=0:1:100
 thr = prctile(lam,iprc);
th_list= [th_list;thr]; 
end

meanl=[];meanR=[];
for ith=1:length(th_list)-1
iidx = (lam>th_list(ith)) & (lam<=th_list(ith+1));
meanl = [meanl;mean(lam(iidx))];
meanR = [meanR;mean(spksGen_hr(iidx))];
end

plot(meanl,meanR);
hold on;
 end


[XX,NN] =ecdf(lam);
hold on;
plotyy(meanl,meanR,NN(NN<=max(meanl)),XX(NN<=max(meanl)));

