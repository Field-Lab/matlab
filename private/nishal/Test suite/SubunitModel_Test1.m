% Null space simulation suite ! 

%% Generate RFs - subunit model
nSubunits = 4;
Filtdim1 = 6;
Filtdim2 = 6;
Filtlen = 30;

subunits=cell(nSubunits,1);

k=zeros(Filtdim1,Filtdim2,1,Filtlen);
k(320-318,160-158,:,:)=1;
k(323-318,162-158,:,:)=1;
k(323-318,163-158,:,:)=1;
subunit_scale=1%1;
subunits{1}=k*subunit_scale;

k=zeros(Filtdim1,Filtdim2,1,Filtlen);
k(320-318,159-158,:,:)=1;
subunit_scale=1%1.2;
subunits{2}=k*subunit_scale;

k=zeros(Filtdim1,Filtdim2,1,Filtlen);
k(322-318,161-158,:,:)=1;
subunit_scale=1%0.5;
subunits{3}=k*subunit_scale;

k=zeros(Filtdim1,Filtdim2,1,Filtlen);
k(322-318,160-158,:,:)=1;
k(319-318,159-158,:,:)=1;
subunit_scale=1%1.5;
subunits{4}=k*subunit_scale;


figure;
for isubunit=1:nSubunits
subplot(2,2,isubunit);
imagesc(subunits{isubunit}(:,:,1,4));
colormap gray
title(sprintf('Subunit: %d',isubunit));
end

figure;
mask=double((subunits{1}(:,:,1,4)+subunits{2}(:,:,1,4)+subunits{3}(:,:,1,4)+subunits{4}(:,:,1,4))~=0);
imagesc(mask);
% Temporal properties?

scale_one=1;
scale_two=0.25;
tau_one=4;
tau_two=10;
n_filters=6;
t=[0:29];
tf = scale_one*((t/tau_one).^n_filters).*exp(-n_filters*(t/tau_one -1)) - scale_two*((t/tau_two).^n_filters).*exp(-n_filters*(t/tau_two -1));
figure;
plot(tf)

tf2=zeros(1,1,1,Filtlen);
tf2(1,1,1,:)=tf;
tf=tf2;
clear tf2
tfRep=repmat(tf,[Filtdim1,Filtdim2,1,1]);
title('Temporal Filter');

subunits{1}=subunits{1}.*tfRep;
subunits{2}=subunits{2}.*tfRep;
subunits{3}=subunits{3}.*tfRep;
subunits{4}=subunits{4}.*tfRep;


% sub-unit weights
subunitWeights=zeros(nSubunits,1);
subunitWeights(1)=1%1;
subunitWeights(2)=1%1.2;
subunitWeights(3)=1%0.7;
subunitWeights(4)=1%1.5;

% sub-unit non-linearity

f= @(x) double(x>0).*(1.6*x);
%f=@(x) exp(0.6*x);
% Ganglion cell non-linearity
%N= @(x) exp(0.15*(x));
N=@(x) double(x>0)*(0.02).*(3.4*x).^2;
%N = @(x) 15./(1+exp(-1.5*(x-5)));


figure;
for itime=1:30
    itime
for isubunit=1:nSubunits
subplot(2,2,isubunit);
imagesc(subunits{isubunit}(:,:,1,itime));
caxis([min(subunits{isubunit}(:)),max(subunits{isubunit}(:))]);
colormap gray
colorbar
%title(sprintf('Subunit: %d',isubunit));

end
%pause
end

%% Generate white noise movie
movieLen=120*30*60;
mov=zeros(Filtdim1,Filtdim2,movieLen);
movie_idx=2; 
if(movie_idx==1)
mov(319-310,158-150,:)=0.5; % Stimulate a sub-unit
mov(320-310,160-150,:)=0.5;
mov=mov+rand(Filtdim1,Filtdim2,movieLen)*0.3;
end

if(movie_idx==2)
mov=double(rand(Filtdim1,Filtdim2,movieLen)>0.5)-0.5;
end

mov(:,:,1:30)=0;

figure;
for itime=40:50
   
imagesc(mov(:,:,itime));
colormap gray
colorbar
axis image
caxis([-0.5,0.5]);
pause(0.01)
end

%% Generate responses
% Calculate filter output for each sub-unit for each frame and calculate
% number of spikes for each frame-bin (binned response) .. So that would be
% used for STA calculation ? 



mov2=zeros(Filtdim1 ,Filtdim2,movieLen+Filtlen-1);
mov2(:,:,Filtlen:movieLen+Filtlen-1)=mov; % Append zeros before the movie
nTrials=100;
SubUnit_Response_test_movie_script

 %% Calculate STA


% My own STA code 
STA=zeros(Filtdim1,Filtdim2,Filtlen);
useTrial=1;
for iframe=30:movieLen
STA=STA+mov(:,:,iframe:-1:iframe-Filtlen+1)*binnedResponses(iframe(:,useTrial));
end
STA=STA/sum(binnedResponses(:,useTrial));


figure
 for itime=[1:Filtlen]
 imagesc(squeeze(STA(:,:,itime)));colormap gray
 caxis([min(STA(:)),max(STA(:))]);
 colorbar
 pause
 end
 
 
 %% Experiment 1
 
% Calculate null space stimulus


movieLen=120*30*60;
null_compute_subUnit_test
 
% Generate responses to null movie


% Calculate filter output for each sub-unit for each frame and calculate
% number of spikes for each frame-bin (binned response) .. So that would be
% used for STA calculation ? 


movie_new_len=size(mov_new2,3);
mov2=zeros(Filtdim1 ,Filtdim2,movie_new_len+Filtlen-1);
mov2(:,:,Filtlen:movie_new_len+Filtlen-1)=mov_new2; % Append zeros before the movie
nTrials=1;
SubUnit_Response_test_movie_script
binnedResponseNull=binnedResponses;
psth_null=psth_resp;
time_log_null=timeLog;

movie_new_len=size(mov_orig2,3);
mov2=zeros(Filtdim1 ,Filtdim2,movie_new_len+Filtlen-1);
mov2(:,:,Filtlen:movie_new_len+Filtlen-1)=mov_orig2; % Append zeros before the movie
nTrials=1;
SubUnit_Response_test_movie_script
binnedResponseOrig=binnedResponses;
psth_orig=psth_resp;
time_log_orig = timeLog;

[x1,y1]=plotSpikeRaster(binnedResponseNull'>0,'PlotType','vertline');
[x2,y2]=plotSpikeRaster(binnedResponseOrig'>0,'PlotType','vertline');

figure;
subplot(2,1,1);
plot(x1,y1,'r');
hold on
plot(x2,y2+max(y2),'k');
xlim([0 max(time_log_orig)]);
ylim([0,2*max(y2)]);

subplot(2,1,2);
plot(time_log_null,psth_null,'k');
hold on
plot(time_log_orig,psth_orig,'r');
xlim([0,max(time_log_null)]);
legend('Null','Original');


% Re-STA
binnedResponses=binnedResponseOrig;
reSTC_SubUnit
reSTAOrig=reSTA;
reSTCOrig=reSTC;

binnedResponses=binnedResponseNull;
reSTC_SubUnit
reSTANull=reSTA;
reSTCNull=reSTC;

%
[V,D]=eigs(reSTCOrig,reSTCNull,10,'lm');
figure;
plot(diag(abs(D)),'*');
title('Eigen Values');

uSq=cell(size(V,2),1);
isel=1;
uSq{isel}=reshape(V(:,isel),[Filtdim1,Filtdim2,Filtlen]).*repmat(mask,[1,1,Filtlen]);
figure; 
for itime=1:30
imagesc(squeeze(uSq{isel}(:,:,itime)));
 colormap gray
 caxis([min(uSq{isel}(:)),max(uSq{isel}(:))])
 hold on
pause(1)
end



%% Experiment 2
 
% Calculate null space stimulus


movieLen=120*15;
null_compute_subUnit_test
 
% Generate responses to null movie


% Calculate filter output for each sub-unit for each frame and calculate
% number of spikes for each frame-bin (binned response) .. So that would be
% used for STA calculation ? 


movie_new_len=size(mov_new2,3);
mov2=zeros(Filtdim1 ,Filtdim2,movie_new_len+Filtlen-1);
mov2(:,:,Filtlen:movie_new_len+Filtlen-1)=mov_new2; % Append zeros before the movie
nTrials=100;
SubUnit_Response_test_movie_script
binnedResponseNull=binnedResponses;
psth_null=psth_resp;
time_log_null=timeLog;

movie_new_len=size(mov_orig2,3);
mov2=zeros(Filtdim1 ,Filtdim2,movie_new_len+Filtlen-1);
mov2(:,:,Filtlen:movie_new_len+Filtlen-1)=mov_orig2; % Append zeros before the movie
nTrials=100;
SubUnit_Response_test_movie_script
binnedResponseOrig=binnedResponses;
psth_orig=psth_resp;
time_log_orig = timeLog;

[x1,y1]=plotSpikeRaster(binnedResponseNull'>0,'PlotType','vertline');
[x2,y2]=plotSpikeRaster(binnedResponseOrig'>0,'PlotType','vertline');

figure;
subplot(2,1,1);
plot(x1,y1,'r');
hold on
plot(x2,y2+max(y2),'k');
xlim([0 max(time_log_orig)]);
ylim([0,2*max(y2)]);

subplot(2,1,2);
plot(time_log_null,psth_null,'k');
hold on
plot(time_log_orig,psth_orig,'r');
xlim([0,max(time_log_null)]);
legend('Null','Original');


% Re-STA
binnedResponses=binnedResponseOrig;
reSTC_SubUnit
reSTAOrig=reSTA;
reSTCOrig=reSTC;

binnedResponses=binnedResponseNull;
reSTC_SubUnit
reSTANull=reSTA;
reSTCNull=reSTC;

%
[V,D]=eigs(reSTCOrig,reSTCNull,10,'lm');

figure;
plot(diag(abs(D)),'*');
title('Eigen Values');


uSq=cell(size(V,2),1);
isel=1;
uSq{isel}=reshape(V(:,isel),[Filtdim1,Filtdim2,Filtlen]).*repmat(mask,[1,1,Filtlen]);
figure; 
for itime=1:30
imagesc(squeeze(uSq{isel}(:,:,itime)));
 colormap gray
 caxis([min(uSq{isel}(:)),max(uSq{isel}(:))])
 hold on
pause(1)
end



%%
% % addpath('../code_stc');
% % 
% %  mov_stc=zeros(size(mov,3),Filtdim1*Filtdim2);
% %  for itime=1:size(mov,3)
% %  xx=squeeze(mov(:,:,itime));
% %  mov_stc(itime,:)=xx(:)';
% %  end
% %  
% %  
% %  [STA,STC] = simpleSTC(mov_stc, binnedResponses, 50);
% %  STAr=reshape(STA,[30,Filtdim1,Filtdim2]);
% %  
% %  figure
% %  for itime=[1:50]
% %  imagesc(squeeze(STAr(itime,:,:)));colormap gray
% %  caxis([min(STAr(:)),max(STAr(:))]);
% %  colorbar
% %  pause
% %  end
% % %  figure;
% % %  [u,s,v]=svds(STC,5);
% % %  plot(diag(s),'*')
% % %  
% % %  uSq=reshape(u(:,2),[30,Filtdim1,Filtdim2]);
% % %  
% % %  figure
% % %  imagesc(squeeze(uSq(4,:,:)));
% % %  colormap gray

% Take responses

% calculate re-STA/STC .. 