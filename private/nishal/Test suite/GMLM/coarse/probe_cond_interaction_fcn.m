
function sd = probe_cond_interaction_fcn(WN_datafile,movie_xml,stim_length,cellID,user_STA_depth,ASM_link,nSU,suP,u1,u2,upossible,sus)
%% Load data


datarun=load_data(WN_datafile)
datarun=load_params(datarun)
datarun=load_sta(datarun);
datarun=load_neurons(datarun);

%% Find stimulus and response
% cellID=1531;

% user_STA_depth=30;
extract_movie_response2;
%% Load fits
data1 = load(ASM_link);
fitGMLM_log = data1.fitGMLM_log;

mask2=data1.mask;
mask=mask2(:);
fitGMLM = fitGMLM_log{nSU};

%% 
 mask=logical(mask);
K=zeros(sum(mask(:)),nSU);
B = ones(nSU,1);
for isu=1:nSU
K(:,isu) = fitGMLM.Linear.filter{isu};
B(isu) = exp(fitGMLM.Linear.bias{isu});
end


figure;
[B_use,h]= plotSU_withcells(K,logical(mask2),logical(mask2(:)),B);
suptitle('All SU');

%sus = input('Which sus?');
figure;
isu=sus(1);
ss = 0*mask2;ss(logical(mask2))=K(:,isu);ss=ss';
plot_SU_with_idx2(ss,logical(mask2));
title(sprintf('SU: %d',isu));

figure
jsu = sus(2);ss = 0*mask2;ss(logical(mask2))=K(:,jsu);ss=ss';
plot_SU_with_idx2(ss,logical(mask2));

title(sprintf('SU: %d',jsu));

Y_C = spksGen';mov = maskedMovdd;

%% Select movie segments where probe movie is strong

uProbe = zeros(sum(mask(:)),1);
cell_choose_num=[1];jcell=1;total_mask_log = mask(:);

%suP = input('Testing which pixels in SU?')
q = -1*total_mask_log(:,cell_choose_num(jcell))*(sum(total_mask_log(:,cell_choose_num(jcell)))-length(suP)); uProbe = uProbe + q(mask);
uProbe(suP)=length(suP);
uProbe = uProbe/norm(uProbe);
 zz=zeros(80,40);zz(mask) = uProbe;figure;imagesc(zz');colormap gray;axis image
 
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

 figure;
[B_use,h]= plotSU_withcells(sta,logical(mask2),logical(mask2(:)),ones(nSU,1));
axis image
title('STA of stimuli chosen for analysis')




%% Do tests
%  u1 = input('Test SU')%[8,9,18,19,30];
%  u2 = input('Other pixels')%[10,20,32,31];
%  upossible=input('Possible other pixels in SU')

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
 sd=testSU_interaction(usu1,usu2,mask(:),mov,Y_C,jcell,probe,thperc,pts,type,'k')
 hold on;
 
%  
%  subplot(1,2,2);
%  thperc=25;pts=7; type='leq';
%  usu1 = u1;usu2=[];
%  testSU_interaction(usu1,usu2,mask(:),mov,Y_C,jcell,probe,thperc,pts,type,'b')
%  hold on;
% 
%  u1perm = perms(u1);u2perm=perms(u2); nperms=size(u2perm,1);npix = length(u1);
%  for iperm=1:nperms
%      
%  usu1 = u1;
%  usu2 = u2perm(iperm,1:(npix));
%   testSU_interaction(usu1,usu2,mask(:),mov,Y_C,jcell,probe,thperc,pts,type,'r')
%  hold on;
%  end
%  
%  usu1 = u1;
%  usu2 = upossible;
%  sd = testSU_interaction(usu1,usu2,mask(:),mov,Y_C,jcell,probe,thperc,pts,type,'k');
%  hold on;
% 


 
 %% 
end