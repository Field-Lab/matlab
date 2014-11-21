data = [2007-8-24-4	10	2	4	8	48.2458	0.4746	48.2118	0.3992	48.2525	0.8105
2007-8-24-4	11	1	4	16	48.1782	0.6345	48.2217	0.3242	48.6505	1.2637
2007-8-24-4	2	9	8	8	96.2906	0.7759	97.0162	1.1235	96.0080	1.4017
2007-8-24-4	6	3	8	8	96.3778	0.8529	96.6834	1.6250	95.8793	1.7374
2007-03-27-1	12	14	8	8	94.2439	1.4895	92.6371	1.6115	95.7292	1.5379
2007-03-27-1	15	3	8	8	93.3637	1.6917	92.8020	1.8534	95.1667	1.6961
2007-8-24-4	7	3	8	16	96.2764	0.8370	96.8241	1.0359	96.0624	1.9367
2007-03-27-1	18	1	8	16	93.4745	1.4377	90.9224	1.4474	93.5113	1.7181
2007-8-24-4	4	1	16	8	192.5311	5.2138	205.3454	19.8429	190.0692	31.1278
2007-03-27-1	16	4	16	8	230.7718	12.7627	218.8598	21.4729	227.6223	26.8721
2007-8-24-4	5	4	16	16	192.3259	3.5522	195.8779	5.8958	188.3950	7.8328
2007-8-24-4	3	23	16	16	191.3860	2.8556	195.3012	5.5379	189.7446	6.1453
2007-03-27-1	13	1	16	16	234.4626	5.5423	224.4872	7.9576	233.5569	8.4548
2007-03-27-1	19	1	16	16	229.9523	5.9215	221.9852	12.0049	234.4199	15.7172];

meanOnP = data(:,6);
stdOnP = data(:,7);

meanOffP = data(:,8);
stdOffP = data(:,9);

meanOnM = data(:,10);
stdOnM = data(:,11);
 filter = data(:,5);
 delta = data(1:8,4)*12;
 date = zeros(8,1);
%  date = repmat('2007-08-24-4',8,1);
 for i = 1:sum(data(:,1) == 1976)
     ind = find(data(:,1) == 1976)
%  date(ind(i),:) = '2007-03-27-1';
  date(ind(i),:) = 1;

 end
 
stdOnP = stdOnP(1:8,:)./delta;
stdOffP = stdOffP(1:8,:)./delta;
stdOnM = stdOnM(1:8,:)./delta;
date = date(1:8,:);

%  figure; 
%  plot(stdOnP, stdOnM,'ko', 'MarkerFaceColor', 'k')
%  hold on
 figure
semilogy(meanOnP(date ==0), stdOnP(date ==0), '*r')
hold on
hs = semilogy(meanOffP(date ==0),stdOffP(date ==0), '*b');
semilogy(meanOnM(date ==0), stdOnM(date ==0), '*k')
semilogy(meanOnP(date ==1), stdOnP(date ==1), 'or')
semilogy(meanOffP(date ==1),stdOffP(date ==1), 'ob');
semilogy(meanOnM(date ==1), stdOnM(date ==1), 'ok')


% figure;
% plot(delta, stdOnP,'.r')
% hold on
% plot(delta, stdOffP,'b.')
% plot(delta, stdOnM, 'g.');
% xlim([3 10])
gap = [0.1 0.1];
marg_h = [.1 .1];
marg_w = [.1 .1];
axes_lim = [0 0.03];
hFig = figure;
set(hFig, 'Position', [680 680 600*1.4 435*1.4])
% ha = tight_subplot_c(2, 3, gap, marg_h, marg_w)

subplot_c(2,3,4)
plot(stdOnP(delta==48), stdOffP(delta==48), 'bo', 'MarkerFaceColor','b')
hold on
plot(stdOnP(delta==96), stdOffP(delta==96), 'ro', 'MarkerFaceColor','r')
axis equal
axis([axes_lim, axes_lim])
hold on
plot(axes_lim,axes_lim, '--k', 'LineWidth',2)
xlabel('ON Parasol STD/Speed')
ylabel('OFF Parasol STD/Speed')

s = sprintf('%c', char(176));
legend(['Speed of 10.7' s '/sec'],['Speed of 21.3' s '/sec'])
subplot_c(2,3,5)
plot(stdOnP(delta==48), stdOnM(delta==48), 'bo', 'MarkerFaceColor','b')
hold on
plot(stdOnP(delta==96), stdOnM(delta==96), 'ro', 'MarkerFaceColor','r')

axis equal
axis([axes_lim, axes_lim])
hold on
plot(axes_lim,axes_lim, '--k', 'LineWidth',2)
xlabel('ON Parasol STD/Speed')
ylabel('ON Midget STD/Speed')


subplot_c(2,3,6)
plot(stdOffP(delta==48), stdOnM(delta==48), 'bo', 'MarkerFaceColor','b')
hold on
plot(stdOffP(delta==96), stdOnM(delta==96), 'ro', 'MarkerFaceColor','r')

axis equal
axis([axes_lim, axes_lim])

hold on
plot(axes_lim,axes_lim, '--k', 'LineWidth',2)
xlabel('OFF Parasol STD/Speed')
ylabel('ON Midget STD/Speed')



subplot_c(2,3,[1 3])
stdOnP_8 = stdOnP(1:8,:);
stdOffP_8 = stdOffP(1:8,:);
stdOnM_8 = stdOnM(1:8,:);
date_8 = date(1:8,:);
% plot(ones(size(stdOnP_8(date_8 ==0))), stdOnP_8(date_8 ==0), 'o', 'MarkerFaceColor','b')
% hold on
% plot(2*ones(size(stdOnP_8(date_8 ==0))), stdOffP_8(date_8 ==0), 'o', 'MarkerFaceColor','b')
% plot(3*ones(size(stdOnP_8(date_8 ==0))), stdOnM_8(date_8 ==0), 'o', 'MarkerFaceColor','b')
title('Low Speed Variablity')
xlabel('OFF Parasol STD/Speed')
ylabel('ON Midget STD/Speed')

plot(ones(size(stdOnP_8(delta==48))), stdOnP_8(delta==48), 'bo', 'MarkerFaceColor','b')
hold on
plot(2*ones(size(stdOnP_8(delta==48))), stdOffP_8(delta==48), 'bo', 'MarkerFaceColor','b')
plot(3*ones(size(stdOnP_8(delta==48))), stdOnM_8(delta==48), 'bo', 'MarkerFaceColor','b')
plot(ones(size(stdOnP_8(delta==96))), stdOnP_8(delta==96), 'ro', 'MarkerFaceColor','r')
hold on
plot(2*ones(size(stdOnP_8(delta==96))), stdOffP_8(delta==96), 'ro', 'MarkerFaceColor','r')
plot(3*ones(size(stdOnP_8(delta==96))), stdOnM_8(delta==96), 'ro', 'MarkerFaceColor','r')

% axis equal
xlim([0.5 3.5])
set(gca, 'xtick', [1,2,3]);
set(gca, 'xticklabel',{'ON Parasol', 'OFF Parasol', 'ON Midget'})

ylabel('STD/Speed')
set(gcf,'color','w');
title('Variability in Speed Estimate in Different Cell Types');
groups={[1.5,3]};
H=sigstar(groups,[0.009]);

% GROUPS - a cell array defining the pairs of groups to compare. Groups defined 
%          either as pairs of scalars indicating locations along the X axis or as 
%          strings corresponding to X-tick labels. Groups can be a mixture of both 
%          definition types.
% STATS -  a vector of p-values the same length as GROUPS. If empty or missing it's 
%          assumed to be a vector of 0.05s the same length as GROUPS. Nans are treated
%          as indicating non-significance.
% NSORT -  optional, 0 by default. If 1, then significance markers are plotted in 
%          the order found in GROUPS. If 0, then they're sorted by the length of the 
%          bar.

% hsAnno = get(hs, 'Annotation');
% hsLegend = get(hsAnno, 'LegendInformation');
% set(hsLegend, 'IconDisplayStyle', 'off');

% 0 margin, 0.02 (normalized) spacing
% spaceplots(gcf, [0.1 0.1 0.1 0.1], [.03 .03]);

