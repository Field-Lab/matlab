dataSets = masterList();

dataSets = {'2011-12-13-2/data000-0'};

%%
for i=1:length(dataSets)
out_all{i} = batchRetinaSummary(dataSets{i},'on midget','nonlin');
end

%%
out = out_all{1};

% good ones: 3, 6, 8

% plot the subunit nonlinearity
clf
set(gcf,'Position',[630 0 810 600]);
inds = sum(out.nonlin.f.out) ~= 0;
sum(inds)
x = out.nonlin.f.x;
y = mean(out.nonlin.f.out(:,inds),2);
ci = std(out.nonlin.f.out(:,inds),[],2);
subplot(2,2,1);
xlim([-2 2]);
ylim([-1 2.3]);
drawHorzLine(0);
drawVertLine(0);
hold on
plot([-2 2],[-2 2],'Color',[0.5 0.5 0.5],'LineStyle',':');
plot([-2 2],[2 -2],'Color',[0.5 0.5 0.5],'LineStyle',':');
plotErrorBars(x,y,ci,[0.6 0.6 0.6],[0.2 0.2 0.2],0)
xlim([-2 2]);
ylim([-1 2.3]);
set(gca,'XTick',[]);
set(gca,'YTick',[]);
set(gca,'FontSize',20,'FontAngle','italic');
xlabel('Subunit input');
ylabel('Subunit output');
title('Subunit nonlinearity');
box off


% plot the output nonlinearity
inds = sum(out.nonlin.g.out) ~= 0;
x = out.nonlin.g.x;
y = mean(out.nonlin.g.out(:,inds),2);
ci = std(out.nonlin.g.out(:,inds),[],2);
subplot(2,2,3);
xlim([-0.5 7]);
ylim([-0.5 2]);
drawHorzLine(0);
drawVertLine(0);
plotErrorBars(x,y,ci,[0.6 0.6 0.6],[0.2 0.2 0.2],0)
xlim([-0.5 7]);
ylim([-0.5 2]);
set(gca,'XTick',[]);
set(gca,'YTick',[]);
set(gca,'FontSize',20,'FontAngle','italic');
foo = xlabel('Input to spiking nonlinearity');
ylabel('Firing rate');
title('Spiking nonlinearity');
box off


% plot r2 comparison
inds = find(out.nonlin.r2.ln ~= 0);
subplot(2,2,2);
plot([0 0.4],[0 0.4],'Color',[0.3 0.3 0.3]);
hold on
for i=1:length(inds)
  circle(out.nonlin.r2.ln(inds(i)),out.nonlin.r2.sub(inds(i)),...
  0.02,[0.3 0.3 0.3],[0 0 0],1);
end
set(gca,'TickDir','out');
axis equal
set(gca,'XTick',[0 0.4]);
set(gca,'YTick',[0 0.4]);
set(gca,'FontSize',20,'FontAngle','italic');
xlim([0 0.4]);
ylim([0 0.4]);
xlabel('LN fit (R2)');
ylabel('Subunit fit (R2)');
title('Model performance (all stimuli)');
box off


% plot r2 comparison with stim selection
inds = find(out.nonlin.r2select.ln ~= 0);
subplot(2,2,4);
plot([0 0.4],[0 0.4],'Color',[0.3 0.3 0.3]);
hold on
for i=1:length(inds)
  circle(out.nonlin.r2select.ln(inds(i)),out.nonlin.r2select.sub(inds(i)),...
  0.02,[0.3 0.3 0.3],[0 0 0],1);
end
set(gca,'TickDir','out');
axis equal
set(gca,'XTick',[0 0.4]);
set(gca,'YTick',[0 0.4]);
set(gca,'FontSize',20,'FontAngle','italic');
xlim([0 0.4]);
ylim([0 0.4]);
xlabel('LN fit (R2)');
ylabel('Subunit fit (R2)');
title('Model performance (selected stimuli)');
box off