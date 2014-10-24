%%
% load the datarun files
load /Users/freeman/Dropbox/Projects/retinaSubunits/data/2011-06-30-0/data004/repeatrun.mat

%%
% get the cell index for the repeat run
ind = find(datarun{1}.cell_ids==dat.rgcId);
repeat_num = datarun{1}.mapping.map_to_data004{ind};
repeat_ind = find(datarun{2}.cell_ids==repeat_num);
%%
% create a PSTH
if ~isempty(repeat_ind)
  REFRESH = 15;
  bin = 1/REFRESH;
  STIMDUR = 148*bin;
  nReps = 100;
  psth = zeros(148,nReps);
  for irep = 1:nReps
    inds = 0:bin:STIMDUR;
    t = histc(datarun{2}.trials.rasters{repeat_ind}{irep},inds);
    t = t(1:end-1);
    psth(:,irep) = t;
  end
  meanResp = mean(psth,2);
end
%%
% plot the rasters
colormap(gray)
subplot(3,1,1);
imagesc(inds(1:end-1),[1:100],clip(2-psth',0,2));
box off
set(gca,'YDir','normal');
ylabel('Repeats');
xlim([0 10]);

subplot(3,1,2);
plot(inds(1:end-1),meanResp*(1/bin),'k');
box off;
xlabel('Time (Seconds)');
ylabel('Impulses per second (ips');
xlim([0 10]);
%%
% compute a quick sta
shiftAmt = -1;
coneInputs = datarun{2}.cones.cone_inputs(:,dat.coneIds);
mnStim = mean(coneInputs); % get mean stimulus
coneInputs = bsxfun(@minus,coneInputs,mnStim); % subtract it off
% compute sta as a weighted average of stimuli
sta = coneInputs'*circshift(meanResp,[shiftAmt 0]);
subplot(3,2,5);
plot(sta/norm(sta),'b');
hold on
plot(dat.staFit/norm(dat.staFit),'r');
hold off
xlabel('Cone index');
ylabel('STA weights');
ylim([-0.5 0.5]);
subplot(3,2,6);
plot(sta/norm(sta),dat.staFit/norm(dat.staFit),'.');
xlabel('STA, repeats');
ylabel('STA, white noise');
mn = min(sta(:)/norm(sta));
mx = max(sta(:)/norm(sta));
axis equal
xlim([mn mx]);
ylim([mn mx]);
%% GET MODEL FIT - IGNORE FOR NOW
% evaluate model fit
%tmp.X_ct = zscore(coneInputs)';
%tmp.R_t = zeros(1,size(tmp.X_ct,2));
%out = evalFit(test,res.fit_LN);
%subplot(3,1,2);
%hold on
%plot(inds(1:end-1),out.Z_t*30,'b');
%hold off