function [ThresholdCurrent,V,G,Gd,D2,D2d]=NS_FindFieldPropertiesForThreshold(DataPath,ClusterFileName,Pattern,Movies,Channels,x,y,z,direction);

NS_GlobalConstants=NS_GenerateGlobalConstants(61);

[Amplitudes0,Curves0]=NS_ReadCurvesFromClusterFileNew(DataPath,ClusterFileName,Channels,Pattern,Movies,NS_GlobalConstants);
a=abs(Amplitudes0(1,:));
s1=NS_CorrectCurve(Curves0(1,:));
ThresholdCurrent=NS_FindThresholdCurrentBimodalFit(a,s1,50);
Amplitudes1=Amplitudes0(:,1)./Amplitudes0(1,1)*ThresholdCurrent;
[GX,GY,GZ]=NS_EstimatePotentialGradientForPattern(Amplitudes1,Channels,x,y,z);
[D2X,D2Y,D2Z]=NS_EstimatePotentialSecondDevForPattern(Amplitudes1,Channels,x,y,z);
[V]=NS_EstimatePotential(Amplitudes1,Channels,x,y,z);

G=sqrt(GX^2+GY^2+GZ^2);
D2=sqrt(D2X^2+D2Y^2+D2Z^2);

Gd=[GX GY GZ]*direction';
D2d=[D2X D2Y D2Z]*direction';