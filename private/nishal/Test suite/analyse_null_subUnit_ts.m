% Generate responses to null movie


% Calculate filter output for each sub-unit for each frame and calculate
% number of spikes for each frame-bin (binned response) .. So that would be
% used for STA calculation ? 

% Response to null movie
movie_new_len=size(mov_new2,3);
mov2=zeros(Filtdim1 ,Filtdim2,movie_new_len+Filtlen-1);
mov2(:,:,Filtlen:movie_new_len+Filtlen-1)=mov_new2; % Append zeros before the movie
SubUnit_Response_test_movie_script
binnedResponseNull=binnedResponses;
psth_null=psth_resp;
time_log_null=timeLog;
cell_resp_orig=cell_resp;

% Response to original movie
movie_new_len=size(mov_orig2,3);
mov2=zeros(Filtdim1 ,Filtdim2,movie_new_len+Filtlen-1);
mov2(:,:,Filtlen:movie_new_len+Filtlen-1)=mov_orig2; % Append zeros before the movie
SubUnit_Response_test_movie_script
binnedResponseOrig=binnedResponses;
psth_orig=psth_resp;
time_log_orig = timeLog;
cell_resp_null=cell_resp;

figure;


plot(max(cell_resp_orig(uint16(cell_resp_len/2-100):uint16(cell_resp_len/2+100),:)'),'b');
hold on
plot(max(cell_resp_null(uint16(cell_resp_len/2-100):uint16(cell_resp_len/2+100),:)'),'b--');
hold on

plot(min(cell_resp_orig(uint16(cell_resp_len/2-100):uint16(cell_resp_len/2+100),:)'),'r');
hold on
plot(min(cell_resp_null(uint16(cell_resp_len/2-100):uint16(cell_resp_len/2+100),:)'),'r--');
hold on

plot(max(abs(cell_resp_orig(uint16(cell_resp_len/2-100):uint16(cell_resp_len/2+100),:))'),'k');
hold on
plot(max(abs(cell_resp_null(uint16(cell_resp_len/2-100):uint16(cell_resp_len/2+100),:))'),'k--');
hold on

legend('max Orig','max Null','min Orig','min Null','absolute max Orig','absolute max Null','Location','best');
title('Max and min inputs to sub-units in Original and Null movies');

figure
[x1,y1]=plotSpikeRaster(binnedResponseNull'>0,'PlotType','vertline');
[x2,y2]=plotSpikeRaster(binnedResponseOrig'>0,'PlotType','vertline');

figure;
subplot(2,1,1);
plot(x1,y1,'k');
hold on
plot(x2,y2+max(y2),'r');
xlim([0 max(time_log_orig)]);
ylim([0,2*max(y2)]);
title('Rasters');
legend('Null','Original');

subplot(2,1,2);
plot(time_log_null,psth_null,'k');
hold on
plot(time_log_orig,psth_orig,'r');
xlim([0,max(time_log_null)]);
legend('Null','Original');
title('PSTH')
% 
% figure;
% scatter(psth_orig',psth_null');
% title('Scatter between Original PSTH and null PSTH');

% figure;
% scatter(binnedResponseOrig,binnedResponseNull);
% title('Scatter between Original and Null response');

%% Re-STA and re STC
binnedResponses=binnedResponseOrig;
reSTC_SubUnit_subtractSTA
reSTAOrig=reSTA;
reSTCOrig=reSTC;

binnedResponses=binnedResponseNull;
reSTC_SubUnit_subtractSTA
reSTANull=reSTA;
reSTCNull=reSTC;

% 
[V,D]=eigs(reSTCNull,reSTCOrig,10,'lm');
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
