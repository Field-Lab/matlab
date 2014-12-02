% data = [2007-8-24-4	10	2	4	8	48.2458	0.4746	48.2118	0.3992	48.2525	0.8105
%bright bars moving right
data = [2007-8-24-4	10	2	4	8	48.2458	0.4746	48.2118	0.3992	48.2525	0.8105
2007-8-24-4	11	1		4	16	48.1782	0.6345	48.2217	0.3242	48.6505	1.2637
2007-8-24-4	 2	9		8	8	96.2906	0.7759	97.0162	1.1235	96.0080	1.4017
2007-8-24-4	 6	3		8	8	96.3778	0.8529	96.6834	1.6250	95.8793	1.7374
2007-03-27-1	 12	14		8	8	117.8055	1.8617	115.7967	2.0145	119.6618	1.9221
2007-03-27-1	 15	3		8	8	116.7047	2.1148	116.0024	2.3161	118.9582	2.12
2007-8-24-4	 7	3		8	16	96.2764	0.8370	96.8241	1.0359	96.0624	1.9367
2007-03-27-1 18	1		8	16	116.8437	1.7976	113.6532	1.8098	116.8896	2.1475
2007-8-24-4	 4	1		16	8	192.5311	5.2138	205.3454	19.8429	190.0692	31.1278
2007-03-27-1 16 	4		16	8	230.7718	12.7627	218.8598	21.4729	227.6223	26.8721
2007-8-24-4	 5 	4		16	16	192.3259	3.5522	195.8779	5.8958	188.3950	7.8328
2007-8-24-4	 3	23		16	16	191.3860	2.8556	195.3012	5.5379	189.7446	6.1453
2007-03-27-1 13	1		16	16	234.4626	5.5423	224.4872	7.9576	233.5569	8.4548
2007-03-27-1 19	1		16	16	229.9523	5.9215	221.9852	12.0049	234.4199	15.7172];
%dark bars moving right
% data = [2007-03-27-1	 12	2	8	8	93.1327	2.5905	93.9007	0.7951	94.3198	1.3409
% 2007-03-27-1	 15	4		8	8	90.8243	3.2959	92.7500	0.7907	93.8439	2.1854
% 2007-03-27-1	 18	4		8	16	91.0135	2.6286	92.2978	0.6092	91.9926	2.5153
% 2007-03-27-1	 16 	1		16	8	184.1164	24.5050	182.9856	5.3567	192.1703	24.8277
% 2007-03-27-1	 13	2		16	16	181.6240	7.7479	183.4271	3.5344	183.6619	10.6257
% 2007-03-27-1	 19	4		16	16	179.9152	10.0789	181.9082	5.0169	181.3347	13.4602
% 2007-8-24-4	 10	4		4	8	48.0354	1.2482	48.3788	0.1816	48.9651	3.5017
% 2007-8-24-4	 11	4		4	16	48.4958	0.4048	48.2571	0.2206	49.5562	4.1184
% 2007-8-24-4	2	7		8	8	96.9170	1.0014	96.6715	0.5256	96.1238	1.7942
% 2007-8-24-4 6	4		8	8	97.3693	1.4487	96.7375	0.4937	95.9965	1.9912
% 2007-8-24-4	7	1		8	16	96.4753	0.7539	97.1396	0.4863	96.6462	2.4552
% 2007-8-24-4	 4	3	16	8	193.8624	8.0861	194.9085	4.2517	193.9645	24.2768
% 2007-8-24-4 5 	2		16	16	193.5828	5.0253	193.7174	2.3587	193.6413	8.5866
% 2007-8-24-4	3	22		16	16	194.4410	4.6275	193.8906	2.1418	198.3383	8.8379];

meanOnP = data(:,6);
stdOnP = data(:,7);

meanOffP = data(:,8);
stdOffP = data(:,9);

meanOnM = data(:,10);
stdOnM = data(:,11);
 filter = data(:,5);
 delta = data(:,4)*12;
 date = zeros(14,1);
%  date = repmat('2007-08-24-4',8,1);
 for i = 1:sum(data(:,1) == 1976)
     ind = find(data(:,1) == 1976)
%  date(ind(i),:) = '2007-03-27-1';
  date(ind(i),:) = 1;

 end
 
stdOnP = stdOnP(:,:)./delta;
stdOffP = stdOffP(:,:)./delta;
stdOnM = stdOnM(:,:)./delta;
date = date(:,:);


%  figure
% semilogy(meanOnP(date ==0), stdOnP(date ==0), '*r')
% hold on
% hs = semilogy(meanOffP(date ==0),stdOffP(date ==0), '*b');
% semilogy(meanOnM(date ==0), stdOnM(date ==0), '*k')
% semilogy(meanOnP(date ==1), stdOnP(date ==1), 'or')
% semilogy(meanOffP(date ==1),stdOffP(date ==1), 'ob');
% semilogy(meanOnM(date ==1), stdOnM(date ==1), 'ok')
% 
% 


axes_lim = [0 0.2];
hFig = figure;
set(hFig, 'Position', [680 680 600*1.5 435*1.5])

subplot_c(2,3,4)
plot(stdOnP(delta==192), stdOffP(delta==192), 'go', 'MarkerFaceColor','g')

hold on
plot(stdOnP(delta==96), stdOffP(delta==96), 'ro', 'MarkerFaceColor','r')
plot(stdOnP(delta==48), stdOffP(delta==48), 'bo', 'MarkerFaceColor','b')

axis equal
axis([axes_lim, axes_lim])
hold on
plot(axes_lim,axes_lim, '--k', 'LineWidth',2)
xlabel('ON Parasol STD/Speed')
ylabel('OFF Parasol STD/Speed')

s = sprintf('%c', char(176));
legend(['Speed of 10.7' s '/sec'],['Speed of 21.3' s '/sec'], ['Speed of 42.6' s '/sec'])
subplot_c(2,3,5)
plot(stdOnP(delta==192), stdOnM(delta==192), 'go', 'MarkerFaceColor','g')

hold on
plot(stdOnP(delta==96), stdOnM(delta==96), 'ro', 'MarkerFaceColor','r')
plot(stdOnP(delta==48), stdOnM(delta==48), 'bo', 'MarkerFaceColor','b')

axis equal
axis([axes_lim, axes_lim])
hold on
plot(axes_lim,axes_lim, '--k', 'LineWidth',2)
xlabel('ON Parasol STD/Speed')
ylabel('ON Midget STD/Speed')


subplot_c(2,3,6)
plot(stdOffP(delta==192), stdOnM(delta==192), 'go', 'MarkerFaceColor','g')
hold on
hold on
plot(stdOffP(delta==96), stdOnM(delta==96), 'ro', 'MarkerFaceColor','r')
plot(stdOffP(delta==48), stdOnM(delta==48), 'bo', 'MarkerFaceColor','b')

axis equal
axis([axes_lim, axes_lim])

hold on
plot(axes_lim,axes_lim, '--k', 'LineWidth',2)
xlabel('OFF Parasol STD/Speed')
ylabel('ON Midget STD/Speed')



subplot_c(2,3,[1 3])
title('Low Speed Variablity')
xlabel('OFF Parasol STD/Speed')
ylabel('ON Midget STD/Speed')
plot(3*ones(size(stdOnP(delta==192))), stdOnM(delta==192), 'go', 'MarkerFaceColor','g')
hold on
plot(ones(size(stdOnP(delta==192))), stdOnP(delta==192), 'go', 'MarkerFaceColor','g')
plot(2*ones(size(stdOnP(delta==192))), stdOffP(delta==192), 'go', 'MarkerFaceColor','g')

plot(ones(size(stdOnP(delta==96))), stdOnP(delta==96), 'ro', 'MarkerFaceColor','r')

hold on
plot(2*ones(size(stdOnP(delta==96))), stdOffP(delta==96), 'ro', 'MarkerFaceColor','r')
plot(3*ones(size(stdOnP(delta==96))), stdOnM(delta==96), 'ro', 'MarkerFaceColor','r')


plot(ones(size(stdOnP(delta==48))), stdOnP(delta==48), 'bo', 'MarkerFaceColor','b')
hold on
plot(2*ones(size(stdOnP(delta==48))), stdOffP(delta==48), 'bo', 'MarkerFaceColor','b')

plot(3*ones(size(stdOnP(delta==48))), stdOnM(delta==48), 'bo', 'MarkerFaceColor','b')


xlim([0.5 3.5])
ylim(axes_lim)
set(gca, 'xtick', [1,2,3]);
set(gca, 'xticklabel',{'ON Parasol', 'OFF Parasol', 'ON Midget'})

ylabel('STD/Speed')
set(gcf,'color','w');
title({'Variability in Speed Estimate in Different Cell Types';'Bright Bar Moving Right'});
% groups={[1.5,3]};
% H=sigstar(groups,[0.009]);


