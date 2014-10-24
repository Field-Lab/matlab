
%% get summaries for lots of data sets
dataSets = masterList();
iCount = 1;
for i=1:length(dataSets.midget.off)
  out{iCount} = batchRetinaSummary(dataSets.midget.off{i},'off midget','counts');
  iCount = iCount + 1;
end

% plot the wassle data
set(gcf,'Position',[0 0 400 400]);
hold on
clrs = cmap_hsv_dark(10);
wassle.ecc = [0.7,2.1; 6,8; 14,15.7];
wassle.subCount = [98 4 0 0 0 0; 38 6 1 0 0 0; 0 6 17 5 0 0];
for iecc = 1:3
  nTotal = sum(wassle.subCount(iecc,:));
  for icount = 3:-1:1
    count = wassle.subCount(iecc,icount);
    if count > 0 || icount < 2
      h = plot(wassle.ecc(iecc,:),repmat(count/nTotal,[1 2]));
      set(h,'LineWidth',20);
      set(h,'Color',clrs(icount,:));
    end
  end
end

% plot the model fits
for i = 1:length(out)
  ecc = getEcc(dataSets.names{i});
  nTotal = nansum(out{i}.counts.subCount(:));
  for icount = 3:-1:1
    count = nansum(out{i}.counts.subCount(:,icount));
    if icount == 1
      foo(i,2) = count/nTotal;
      foo(i,1) = ecc;
    end
    if count > 0
      h = plot(ecc,count/nTotal,'o');
      set(h,'MarkerSize',20);
      set(h,'MarkerFaceColor',clrs(icount,:));
      set(h,'MarkerEdgeColor','w');
    end
  end
end

% fix plots
xlim([-0.3 16.6]);
ylim([-0.1 1.1]);
set(gca,'XTick',[0:5:15]);
set(gca,'YTick',[0:0.5:1]);
set(gca,'FontSize',16);
set(gca,'TickDir','out');
xlabel('Eccentricity (mm)');
ylabel('Probability of occurance')


clf
hold on
for i=1:length(out)
  thisecc = getEcc(dataSets.names{i});
  frac = sum(out{i}.counts.density(:,2))/sum(out{i}.counts.density(:,1));
  foo(i,2) = frac;
  foo(i,1) = thisecc;
  thisecc
  frac
  h = plot(thisecc,frac,'o');
  set(h,'MarkerSize',20);
  set(h,'MarkerFaceColor',[0.5 0.5 0.5]);
end
xlim([4 16]);
ylim([-0.1 1.1]);

b_cone = fitLine(log(ecc.cone.x),ecc.cone.y);
b_bipolar = fitLine(log(ecc.bipolar.x),ecc.bipolar.y);
x_eval = linspace(5,15,100);
eval_cone = b_cone(2) + b_cone(1)*log(x_eval);
eval_bipolar = b_bipolar(2) + b_bipolar(1)*log(x_eval);
plot(x_eval,eval_bipolar./eval_cone);


