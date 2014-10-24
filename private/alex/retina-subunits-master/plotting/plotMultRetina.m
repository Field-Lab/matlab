%% get summaries for lots of data sets
dataSets = masterList();
iCount = 1;
for i=1:length(dataSets.midget.off)
  out{iCount} = batchRetinaSummary(dataSets.midget.off{i},'off midget',{'nonlin','counts'});
  iCount = iCount + 1;
end
for i=1:length(dataSets.midget.on)
  out{iCount} = batchRetinaSummary(dataSets.midget.on{i},'on midget',{'nonlin','counts'});
  iCount = iCount + 1;
end

% for all offs, and then all ons, plot improvement in model performance
% across cells

set(gcf,'Position',[0 0 800 420]);
subplot(1,2,1)
hold on
improvement = [];
for i=1:length(dataSets.midget.off)
	badinds = out{i}.nonlin.r2.ln < 0.05;
	diffVal = out{i}.nonlin.r2.sub(~badinds)./out{i}.nonlin.r2.ln(~badinds);
	improvement(i) = mean(diffVal);
	h = plot(mean(out{i}.nonlin.r2.ln(~badinds)),mean(out{i}.nonlin.r2.sub(~badinds)),'o');
	set(h,'MarkerFaceColor','k');
	set(h,'MarkerEdgeColor','none');
	set(h,'MarkerSize',20);
end
for i=(length(dataSets.midget.off) + 1):length(dataSets.names)
	badinds = out{i}.nonlin.r2.ln < 0.05;
	diffVal = out{i}.nonlin.r2.sub(~badinds)./out{i}.nonlin.r2.ln(~badinds);
	improvement(i) = mean(diffVal);
	h = plot(mean(out{i}.nonlin.r2.ln(:)),mean(out{i}.nonlin.r2.sub(:)),'o');
	set(h,'MarkerFaceColor',[0.9 0.9 0.9]);
	set(h,'MarkerEdgeColor','none');
	set(h,'MarkerSize',20);
end
plot([-1 0.4],[-1 0.4],'k-');
axis square;
xlim([-0.025 0.425]);
ylim([-0.025 0.425]);
set(gca,'XTick',[0 0.4]);
set(gca,'YTick',[0 0.4]);
set(gca,'FontSize',16);
xlabel('LN fit (R2)');
ylabel('Subunit fit (R2)');
set(gca,'TickDir','out');

subplot(1,2,2)
hold on
for i=1:length(dataSets.midget.off)
	badinds = out{i}.nonlin.r2select.ln < 0.05;
	diffVal = out{i}.nonlin.r2select.sub(~badinds)./out{i}.nonlin.r2select.ln(~badinds);
	improvement_select(i) = mean(diffVal);
	h = plot(mean(out{i}.nonlin.r2select.ln(:)),mean(out{i}.nonlin.r2select.sub(:)),'o');
	set(h,'MarkerFaceColor','k');
	set(h,'MarkerEdgeColor','none');
	set(h,'MarkerSize',20);
end
for i=(length(dataSets.midget.off) + 1):length(dataSets.names)
	badinds = out{i}.nonlin.r2select.ln < 0.05;
	diffVal = out{i}.nonlin.r2select.sub(~badinds)./out{i}.nonlin.r2select.ln(~badinds);
	improvement_select(i) = mean(diffVal);
	h = plot(mean(out{i}.nonlin.r2select.ln(:)),mean(out{i}.nonlin.r2select.sub(:)),'o');
	set(h,'MarkerFaceColor',[0.9 0.9 0.9]);
	set(h,'MarkerEdgeColor','none');
	set(h,'MarkerSize',20);
end
plot([-1 0.4],[-1 0.4],'k-');
axis square;
xlim([-0.025 0.425]);
ylim([-0.025 0.425]);
set(gca,'XTick',[0 0.4]);
set(gca,'YTick',[0 0.4]);
set(gca,'FontSize',16);
xlabel('LN fit (R2)');
ylabel('Subunit fit (R2)');
set(gca,'TickDir','out');

fprintf('OFF cells')
mean(improvement(1:length(dataSets.midget.off)))
std(improvement(1:length(dataSets.midget.off)))
mean(improvement_select(1:length(dataSets.midget.off)))
std(improvement_select(1:length(dataSets.midget.off)))

fprintf('ON cells')
mean(improvement(length(dataSets.midget.off)+1:end))
std(improvement(length(dataSets.midget.off)+1:end))
mean(improvement_select(length(dataSets.midget.off)+1:end))
std(improvement_select(length(dataSets.midget.off)+1:end))

