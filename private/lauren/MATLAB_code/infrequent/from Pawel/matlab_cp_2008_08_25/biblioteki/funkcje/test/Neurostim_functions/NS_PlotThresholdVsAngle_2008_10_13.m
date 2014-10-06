function Threshold0=NS_PlotThresholdVsAngle_2008_10_13(DataPath,FileName,CenterChannel,Pattern0,Patterns1,Patterns2,Angles,Movies0,Movies1,Movies2,FigureProperties,NS_GlobalConstants);

[Amplitudes0,Curves0]=NS_ReadCurvesFromClusterFile(DataPath,FileName,CenterChannel,Pattern0,Movies0,NS_GlobalConstants);
[Amplitudes1,Curves1]=NS_ReadCurvesFromClusterFile(DataPath,FileName,CenterChannel,Patterns1,Movies1,NS_GlobalConstants);
[Amplitudes2,Curves2]=NS_ReadCurvesFromClusterFile(DataPath,FileName,CenterChannel,Patterns2,Movies2,NS_GlobalConstants);

%Curves0=NS_CorrectCurve(Curves0);
%Curves0=NS_CorrectCurve(Curves1);
%Curves0=NS_CorrectCurve(Curves1);
s0=NS_CorrectCurve(Curves0(1,:));
Threshold0=NS_FindThresholdCurrent(Amplitudes0,NS_CorrectCurve(Curves0(1,:)),50)

figure(FigureProperties.FigureNumber+100);
AnglesSorted=sort(Angles);
for i=1:6
    s1=NS_CorrectCurve(Curves1(i,:));
    s2=NS_CorrectCurve(Curves2(i,:));
    plot(Amplitudes1,s1,'bd-',Amplitudes2,s2,'rd-');
    hold on;
    Thresholds1(i)=NS_FindThresholdCurrent(Amplitudes1,s1,50);
    Thresholds2(i)=NS_FindThresholdCurrent(Amplitudes2,s2,50);
end

for i=1:6
    angle=AnglesSorted(i);
    a=find(Angles==angle);
    ThresholdsSorted1(i)=Thresholds1(a);
    ThresholdsSorted2(i)=Thresholds2(a);
end

figure(FigureProperties.FigureNumber+200);
plot(Amplitudes0,s0);

AnglesSorted=[AnglesSorted AnglesSorted(1)+360];
ThresholdsSorted1=[ThresholdsSorted1 ThresholdsSorted1(1)];
ThresholdsSorted2=[ThresholdsSorted2 ThresholdsSorted2(1)];
figure(FigureProperties.FigureNumber);
h=plot(AnglesSorted,ThresholdsSorted1/Threshold0,'bd-',AnglesSorted,ThresholdsSorted2/Threshold0,'rd-')
set(h,'LineWidth',2);
axis([0 360 0.6 2]);
h=text(180,1.9,['Threshold on center=' num2str(Threshold0) '\muA']);
set(h,'FontSize',24);
grid on;

h=xlabel('Angle [degrees]');
set(h,'FontSize',24);
h=ylabel('Normalized threshold');
set(h,'FontSize',24);
h=gca;
set(h,'FontSize',24);
