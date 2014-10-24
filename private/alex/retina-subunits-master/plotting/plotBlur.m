%% get summaries for lots of data sets

dataSets = {'2012-04-13-1/data002','2011-10-14-1/data001-0',...
  '2008-08-27-5/data003','2010-03-05-2/data013',...
  '2011-10-25-5/data001-0','2011-12-13-2/data000-0',...
  '2011-07-05-2/data002','2011-10-25-9/data006-0',...
  '2011-06-24-6/data005','2011-06-30-0/data003',...
  '2012-07-26-1/data004','2012-08-09-1/data001',...
  '2012-08-21-0/data001',
  };

%%
dataSets = {'2011-10-14-1/data001-0'};

%%
for i=1:length(dataSets)
  out{i} = batchRetinaSummary(dataSets{i},'off midget','blury');
end

%%
clf
nPlots = (ceil(sqrt(length(dataSets))));
for i=1:length(dataSets)
  subplot(nPlots,nPlots,i);
  x = linspace(0,4,15);
  hsub = hist(cat(2,out{i}.blury.sz.subpairs{:}),x);
  hnonsub = hist(cat(2,out{i}.blury.sz.nonsubpairs{:}),x);
  hold on
  plot(x,hsub/sum(hsub),'r','LineWidth',4);
  plot(x,hnonsub/sum(hnonsub),'k','LineWidth',4);
  [h p] = ttest2(cat(2,out{i}.blury.sz.subpairs{:}),cat(2,out{i}.blury.sz.nonsubpairs{:}))
  %[h p] = kstest2(cat(2,out{i}.blury.sz.subpairs{:}),cat(2,out{i}.blury.sz.nonsubpairs{:}))
  title(sprintf('%s,  p = %.02g',dataSets{i},p),'FontSize',14);
  set(gca,'TickDir','out');
  set(gca,'FontSize',14);
  xlim([0 4]);
  ylim([0 0.75]);
  mn_sub(i) = mean(cat(2,out{i}.blury.sz.subpairs{:}));
  mn_nonsub(i) = mean(cat(2,out{i}.blury.sz.nonsubpairs{:}));
end
%%
i = 1;
rng = [4 8];
dists_sub = cat(1,out{i}.blury.dists.subpairs{:});
normdists_sub = cat(1,out{i}.blury.normdists.subpairs{:});
dists_nonsub = cat(1,out{i}.blury.dists.nonsubpairs{:});
normdists_nonsub = cat(1,out{i}.blury.normdists.nonsubpairs{:});

[~,inds_sub] = isrange(dists_sub,rng(1),rng(2));
mean(normdists_sub(inds_sub))
[~,inds_nonsub] = isrange(dists_nonsub,rng(1),rng(2));
mean(normdists_nonsub(inds_nonsub))
[h p] = ttest2(normdists_sub(inds_sub),normdists_nonsub(inds_nonsub))
%%
clrs = cmap_hsv_dark(10);
clf

hold on
for i = 1:length(out)
  nonsubpairs = cat(1,out{i}.blury.nonsubpairs{:});
  x = ones(1,size(nonsubpairs,1));
  x = x + randn(size(x))*0.1;
  h1 = plot(x,nonsubpairs,'o');
  set(h1,'MarkerSize',10);
  set(h1,'MarkerFaceColor',clrs(1,:));
  set(h1,'MarkerEdgeColor','k');
  
  subpairs = cat(1,out{i}.blury.subpairs{:});
  x = ones(1,size(subpairs,1))*2;
  x = x + randn(size(x))*0.1;
  h2 = plot(x,subpairs,'o');
  set(h2,'MarkerSize',10);
  set(h2,'MarkerFaceColor',clrs(2,:));
  set(h2,'MarkerEdgeColor','k');
end
%%
xlim([0 3]);
xlim([0.5 2.5]);
set(gca,'XTick',[1 2]);
set(gca,'XTickLabel',{'Diff sub','Same sub'});
box off
set(gca,'TickDir','out');
ylim([1.5 5]);
ylabel('Normalized cone pair distance','FontSize',14);
set(gca,'FontSize',14);

%%
xlim([-0.3 16.6]);
ylim([-0.1 1.1]);
xlabel('Eccentricity (mm)');
set(gca,'XTick',[0:5:15]);
set(gca,'YTick',[0:0.25:1]);