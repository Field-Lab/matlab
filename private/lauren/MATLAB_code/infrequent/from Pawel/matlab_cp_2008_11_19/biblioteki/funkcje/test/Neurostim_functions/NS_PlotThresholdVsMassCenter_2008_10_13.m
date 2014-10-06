function Threshold0=NS_PlotThresholdVsMassCenter_2008_10_13(DataPath,FileName,CenterChannel,Pattern0,Patterns,Movies0,Movies,x,y,FigureProperties,NS_GlobalConstants);

[Amplitudes0,Curves0]=NS_ReadCurvesFromClusterFile(DataPath,FileName,CenterChannel,Pattern0,Movies0,NS_GlobalConstants);
[Amplitudes1,Curves1]=NS_ReadCurvesFromClusterFile(DataPath,FileName,CenterChannel,Patterns,Movies,NS_GlobalConstants);

s0=NS_CorrectCurve(Curves0(1,:));
Threshold0=NS_FindThresholdCurrent(Amplitudes0,NS_CorrectCurve(Curves0(1,:)),50);

for i=1:length(Patterns)
    s1=NS_CorrectCurve(Curves1(i,:));
    Thresholds(i)=NS_FindThresholdCurrent(Amplitudes1,s1,50);
    Amplitudes=NS_AmplitudesForPattern(DataPath,Channels,Patterns(i),MovieNumber,NS_GlobalConstants);
    [xs,ys]=NS_EstimateMassCenter(Amplitudes,Channels);    
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
axis([0 360 0.7 1.3]);
h=text(180,1.9,['Threshold on center=' num2str(Threshold0) '\muA']);
set(h,'FontSize',24);
grid on;

h=xlabel('Angle [degrees]');
set(h,'FontSize',24);
h=ylabel('Normalized threshold');
set(h,'FontSize',24);
h=gca;
set(h,'FontSize',24);
