function h=NS_PlotEfficacyVsAngle_2008_10_13(eff0,eff1,eff2,Angles,FigureProperties);

figure(FigureProperties.FigureNumber);

AnglesSorted=sort(Angles);
for i=1:6
    angle=AnglesSorted(i);
    a=find(Angles==angle);
    effSorted1(i)=eff1(a)/eff0*50;
    effSorted2(i)=eff2(a)/eff0*50;
end

AnglesSorted=[AnglesSorted AnglesSorted(1)+360];
effSorted1=[effSorted1 effSorted1(1)];
effSorted2=[effSorted2 effSorted2(1)];
figure(FigureProperties.FigureNumber);
h=plot(AnglesSorted,effSorted1,'bd-',AnglesSorted,effSorted2,'rd-')
set(h,'LineWidth',2);
axis([0 360 0 100]);
%h=text(180,1.9,['Threshold on center=' num2str(eff0) '\muA']);
%set(h,'FontSize',24);
grid on;

h=xlabel('Angle [degrees]');
set(h,'FontSize',24);
h=ylabel('Efficacy [%]');
set(h,'FontSize',24);
h=gca;
set(h,'FontSize',24);
