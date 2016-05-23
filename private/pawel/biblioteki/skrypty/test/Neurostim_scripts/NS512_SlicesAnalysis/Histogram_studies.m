figure(200)
clf
subplot(1,2,1)

p=hist(LatenciesForGivenElectrode,t);
hist(LatenciesForGivenElectrode,t);
CFs=NS512_FitWithMultiGauss(t,p);
%plot(t,p);
hold on;
h=plot(CFs{1})
set(h,'LineWidth',2)
legend off;
h=gca;
set(h,'XLim',[0 600]);
xlabel('Time [ms]');
ylabel('N');
text(450,9,['A=' num2str(CFs{1}.A,'%8.2f')]);
text(450,8,['\tau=' num2str(CFs{1}.tau,'%8.2f')]);
text(450,7,['\sigma=' num2str(CFs{1}.sigma,'%8.2f')]);

subplot(1,2,2);
p=hist(LatenciesForGivenElectrode,[5:10:600]);
hist(LatenciesForGivenElectrode,[5:10:600]);
CFs=NS512_FitWithMultiGauss([1:60],p);
%plot([5:10:600],p);
hold on
FL=CFs{1}
h=plot([5:10:600],FL([1:60]),'r-')
set(h,'LineWidth',2)
h=gca;
set(h,'XLim',[0 600]);
xlabel('Time [ms]');
ylabel('N');
text(450,40,['A=' num2str(CFs{1}.A,'%8.2f')]);
text(450,35,['\tau=' num2str(CFs{1}.tau*10,'%8.2f')]);
text(450,30,['\sigma=' num2str(CFs{1}.sigma*10,'%8.2f')]);
break
FullImageName=['C:\home\Pawel\nauka\analiza\SlicesTTX\2010-09-14-0\analysis_2012_04_10\other_plots\histograms_study.tif'];
h=gcf;
set(h,'PaperUnits','inches');
set(h,'PaperSize',[8 4.5]);
set(h,'PaperPosition',[0 0 8 4.5]); 
print(h, '-dtiff', '-r240', FullImageName);