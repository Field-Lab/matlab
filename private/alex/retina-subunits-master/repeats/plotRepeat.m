function plotRepeat(re,celldat)

% plot repet results for repeat struct re and cell struct celldat
% plot repet results for repeat struct re and cell struct celldat

rng = [15:45];

interval = 5;
refresh_rate = 60;
time_vals = [0:size(re.resp,2)-1]*(interval/refresh_rate);

close all;
set(gcf,'Position',[0 0 700 700]);

subplot(2,1,1);
imagesc(clip(max(re.resp(:))-re.resp,0,max(re.resp(:))));
axis off;
colormap(gray);
set(gca,'FontSize',18);
title(sprintf('Cell #%g, Reliability: %.02g',celldat.rgcId,re.reliability));

subplot(2,1,2);
hold on;
convert_to_firing = 1/(interval/refresh_rate);
plot(time_vals,re.R_t*convert_to_firing,'k','LineWidth',5);
plot(time_vals,re.predic.Z_t*convert_to_firing,'Color',[0.3 0.3 0.6],'LineWidth',5);
plot(time_vals,re.predic_LN.Z_t*convert_to_firing,'Color',[0.2 0.5 0.3],'LineWidth',5);
xlim([0 max(time_vals)]);
set(gca,'XTick',[0 2 4 6 8]);

box off;
set(gca,'FontSize',18,'FontAngle','italic');
set(gca,'TickDir','out');
ylabel('Firing rate (ips)');
xlabel('Seconds');
title(sprintf('Explained variance: %.02g vs %.02g',re.cr_boot/re.reliability,re.cr_LN_boot/re.reliability));

hold on
drawVertLine([min(rng)-1,max(rng)-1]*(interval/refresh_rate));

figure;
set(gcf,'Position',[500 0 700 300]);
subplot(3,1,1);
mx = max(re.resp(:));
imagesc(clip(mx-re.resp(:,rng),0,mx),[0 mx]);
axis off;

subplot(3,1,2);
predic = poissrnd(repmat(mean(re.predic.Z_t,1),size(re.resp,1),1));
imagesc(clip(mx-predic(:,rng),0,mx),[0 mx]);
axis off;

subplot(3,1,3);
predic = poissrnd(repmat(mean(re.predic_LN.Z_t,1),size(re.resp,1),1));
imagesc(clip(mx-predic(:,rng),0,mx),[0 mx]);
axis off;

colormap(gray);