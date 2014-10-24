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


foo = [];
for i=1:length(dataSets.midget.off)
	improvement(i) = fitLine(out{i}.nonlin.r2.ln,out{i}.nonlin.r2.sub,[],0,[]);
	badinds = out{i}.nonlin.r2select.ln < 0;
	improvementSelect(i) = fitLine(out{i}.nonlin.r2select.ln(~badinds),out{i}.nonlin.r2select.sub(~badinds),[],0,[]);
	improvementSing(i) = fitLine(out{i}.nonlin.r2.sing,out{i}.nonlin.r2.sub,[],0,[]);
end

for i=length(dataSets.midget.off):length(out)
	improvementOn(i) = fitLine(out{i}.nonlin.r2.ln,out{i}.nonlin.r2.sub,[],0,[]);
	badinds = out{i}.nonlin.r2select.ln < 0;
	improvementSelectOn(i) = fitLine(out{i}.nonlin.r2select.ln(~badinds),out{i}.nonlin.r2select.sub(~badinds),[],0,[]);
	improvementSingOn(i) = fitLine(out{i}.nonlin.r2.sing,out{i}.nonlin.r2.sub,[],0,[]);
end
improvementOn = improvementOn(length(dataSets.midget.off)+1:end);
improvementSelectOn = improvementSelectOn(length(dataSets.midget.off)+1:end);
improvementSingOn = improvementSingOn(length(dataSets.midget.off)+1:end);


subplot(2,1,1);
x = linspace(0,3,50);
h1 = hist(improvement,x);
h2 = hist(improvementOn,x);
bar(x,h1,'EdgeColor','none');
hold on
bar(x,h2,'EdgeColor','none','FaceColor','r');
box off
ylim([0 7]);
set(gca,'YTick',[0 7]);
set(gca,'XTick',[0 1 2 3]);
xlim([0 3]);
set(gca,'TickDir','out');
drawVertLine(1);

subplot(2,1,2);
x = linspace(0,3,50);
h1 = hist(improvementSelect(improvementSelect<4),x);
h2 = hist(improvementSelectOn,x);
bar(x,h1,'EdgeColor','none');
hold on
bar(x,h2,'EdgeColor','none','FaceColor','r');
box off
ylim([0 4]);
set(gca,'YTick',[0 4]);
set(gca,'XTick',[0 1 2 3]);
xlim([0 3]);
set(gca,'TickDir','out');
drawVertLine(1);




for i=1:length(out)
	keepinds = out{i}.nonlin.r2.sub > 0.05;
	foo = out{i}.counts.subCount(keepinds',:);
	foo2 = out{i}.counts.density(keepinds',1);
	singfrac(i) = sum(foo(:,1) == foo2)/sum(keepinds);
	doubfrac(i) = sum(((foo(:,1) + foo(:,2)*2) == foo2) & ~(foo(:,1) == foo2))/sum(keepinds);
	tripfrac(i) = sum(((foo(:,1) + foo(:,2)*2 + foo(:,3)*3) == foo2) & ~(foo(:,1) == foo2) & ~(foo(:,1) + foo(:,2)*2 == foo2))/sum(keepinds);
end
