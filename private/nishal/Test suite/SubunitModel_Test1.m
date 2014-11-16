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
%N= @(x) double(x>0)*(0.000005).*(30*x).^2;
N = @(x) 15./(1+exp(-1.5*(x-5)));


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
    itime
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
 
 
 %% Calculate null space stimulus


movieLen=120*60*60;
mov_gen2=zeros(Filtdim1,Filtdim2,movieLen);
movie_idx=2;


if(movie_idx==1)
mov_gen2(319-310,158-150,:)=0.5; % Stimulate a sub-unit
mov_gen2(320-310,160-150,:)=0.5;
mov_gen2=mov_gen2+rand(Filtdim1,Filtdim2,movieLen)*0.3;
end

if(movie_idx==2)
mov_gen2=double(rand(Filtdim1,Filtdim2,movieLen)>0.5)-0.5;
end
if(movie_idx==3)
    latency=30;
mov_short=double(rand(Filtdim1,Filtdim2,ceil(movieLen/latency)+2)>0.5)-0.5;
mov_gen2=zeros(Filtdim1,Filtdim2,movieLen);
    for iframe=1:movieLen
        mov_gen2(:,:,iframe)=mov_short(:,:,floor(iframe/latency)+1);
    end
end
if(movie_idx==4)
    latency=20;
mov_short=double(rand(Filtdim1,Filtdim2,ceil(movieLen/latency)+2));
mov_gen2=zeros(Filtdim1,Filtdim2,movieLen);
    for iframe=1:movieLen
        mov_gen2(:,:,iframe)=mov_short(:,:,floor(iframe/latency)+1);
    end
end

mov_gen2(:,:,1:30)=0;

figure;
for itime=40:50
    itime
imagesc(mov_gen2(:,:,itime));
colormap gray
colorbar
axis image
caxis([-0.5,0.5]);
pause(0.01)

end


addpath(genpath('../create_act_2/'));

mov2=zeros(Filtdim1,Filtdim2,movieLen+2*120);
mov2(:,:,121:end-120)=mov_gen2;

stas=cell(1,1);
stas{1}=zeros(Filtdim1,Filtdim2,1,Filtlen);
for itime=1:Filtlen
    itime
stas{1}(:,:,1,itime)=STA(:,:,itime).*mask;% STA from WN run and masked by significant pixels. Assume I have correct mask from significant_stixles calculation%mask*tf(1,1,1,itime);%STA(:,:,itime);
end

 [mov_orig2,mov_new2]=fourier_project_2(stas,mov2);
 
 % See movies
 figure;
 for itime=1000:1020
  subplot(3,2,1);
  imagesc(mov_orig2(:,:,itime));
  caxis([min(mov_orig2(:)),max(mov_orig2(:))]);
  colormap gray
  colorbar
  axis image
 
  subplot(3,2,2);
  imagesc(mov_new2(:,:,itime));
  caxis([min(mov_new2(:)),max(mov_new2(:))]);
  colormap gray
  colorbar
  axis image
  
  subplot(3,2,3);
  imagesc(mov_orig2(:,:,itime).*mask);
  caxis([min(mov_orig2(:)),max(mov_orig2(:))]);
  colormap gray
  colorbar
  axis image
 
  
  subplot(3,2,4);
 
  imagesc(mov_new2(:,:,itime).*mask);
  caxis([min(mov_new2(:)),max(mov_new2(:))]);
  colormap gray
  colorbar
  axis image
 
  subplot(3,2,5);
  imagesc((mov_new2(:,:,itime)-mov_orig2(:,:,itime)));
  colormap gray
  colorbar
  axis image
  
  pause(1/120)
 end
 
 % Movie histograms 
  
  mov_gen2_masked=mov_gen2.*repmat(mask,[1,1,size(mov_gen2,3)]);
  mov_new2_center=mov_new2(:,:,121:end-120);
  mov_new2_masked = mov_new2_center.*repmat(mask,[1,1,size(mov_new2_center,3)]);
  
 [C,X]=hist(mov_gen2_masked(mov_gen2_masked~=0),100);
 [C2,X2]=hist(mov_new2_masked(mov_gen2_masked~=0),100);
 figure;
 subplot(2,1,1);
 plotyy(X,C,X2,C2);
 
 legend('Original','Null space');
 title('Movie pixel histograms in mask');
 
subplot(2,1,2);

difference =  mov_gen2_masked - mov_new2_masked;
hist(difference(mov_gen2_masked~=0),20);
title('Pixel value difference in mask');
 
%% Generate responses to null movie


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


%% Re-STA
step=30;
nTrials=1;
% My own STA code 
reSTA=zeros(Filtdim1,Filtdim2,Filtlen);
icnt=0;
for itrial=1:nTrials
for iframe=59:step:movie_new_len
    if(binnedResponseNull(iframe,itrial)>0)
reSTA=reSTA+mov_new2(:,:,iframe:-1:iframe-Filtlen+1)*binnedResponseNull(iframe,itrial);
icnt=icnt+binnedResponseNull(iframe,itrial);
    end
end
end
reSTA=reSTA/icnt;


figure
 for itime=1:Filtlen
 imagesc(squeeze(reSTA(:,:,itime)));colormap gray
 caxis([min(reSTA(:)),max(reSTA(:))]);
 colorbar
 pause
 end
 
 % Relation between STA and re-STA

 % Histograms
 figure
 [X,C]=hist(STA(:),20);
 [X2,C2]=hist(reSTA(:),20);
 plotyy(C,X,C2,X2);
 title('Histograms of STA and re-STA')
 legend('STA','reSTA')
 
 % Correlation between STA and re-STA
 corrSTA=STA.*reSTA;
 corr_STAr=sum(corrSTA(:))/(norm(STA(:))*norm(reSTA(:))) 
 disp('Obviously, Correlation between STA and reSTA too low as its in null space.. ')
 
 %% Re- STC
 
% My own STA code 



reSTC=zeros(Filtdim1*Filtdim2*Filtlen);
indx=[59:1:movie_new_len];
iTrial=1;
binnedResponses=binnedResponseNull(:,iTrial);



framesValid = indx(binnedResponses(indx)>0);

for iframe=framesValid
    iframe;
  xx=  reshape(mov_new2(:,:,iframe:-1:iframe-Filtlen+1),[Filtdim1*Filtdim2*Filtlen,1]);
  
reSTC=reSTC+ xx*xx'*binnedResponses(iframe);
end
yy=reshape(reSTA,[Filtdim1*Filtdim2*Filtlen,1]);
reSTC=reSTC / (sum(binnedResponses)-1) - yy*yy'* (sum(binnedResponses))/(sum(binnedResponses)-1);

[u,s,v]=svds(reSTC,100);
figure;
plot(diag(s),'*')
%%
 
figure;
for isubunit=1:nSubunits
subplot(2,2,isubunit);
imagesc(subunits{isubunit}(:,:,4));
colormap gray
colorbar
end

uSq=cell(4,1); 
for isub=1:4
uSq{isub}=reshape(u(:,isub),[Filtdim1,Filtdim2,Filtlen]).*repmat(mask,[1,1,Filtlen]);
end


 figure
 
 for itime=4%1:Filtlen
 for isub=1:4
     subplot(2,2,isub);
     imagesc(squeeze(uSq{isub}(:,:,itime)));
 colormap gray
 caxis([min(uSq{isub}(:)),max(uSq{isub}(:))])
 colorbar
 %pause
 end
 end
 
%  figure
%  for itime=1:Filtlen
%      itime
%      for iidx=1:36
%      subplot(6,6,iidx)
%      uSq=reshape(u(:,iidx),[Filtdim1,Filtdim2,Filtlen]).*repmat(mask,[1,1,Filtlen]);
%      imagesc(squeeze(uSq(:,:,itime)));
%      colormap gray
%      %caxis([min(u(:)),max(u(:))]);
%       caxis([min(uSq(:)),max(uSq(:))]);
%      colorbar
%      end
% pause
%  end
%  
 
%uSq=reshape(u(:,72),[Filtdim1,Filtdim2,Filtlen]).*repmat(mask,[1,1,Filtlen]);
cell_resp=zeros(1,nSubunits);

for isubunit=1:nSubunits
   
cell_resp(:,isubunit)=sum(uSq(:).*squeeze(subunits{isubunit}(:)));
end
cell_resp
% 
totalInput=0*cell_resp(:,1);
for isubunit=1:nSubunits
totalInput=totalInput+subunitWeights(isubunit)*f(cell_resp(:,isubunit));
end

totalOutput=N(totalInput)

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