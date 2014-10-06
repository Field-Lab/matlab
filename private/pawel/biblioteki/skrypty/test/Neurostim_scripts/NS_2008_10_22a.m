DataPath='E:\analysis\2008-08-26-0\data008_proba4';
ArtifactDataPath='D:\analysis\2008-08-27-4\data014';
ArtifactSubtraction=0;
FileName='E:\analysis\2008-08-26-0\data008_proba4\clusters008';

DataPath='D:\analysis\2008-08-27-4\data005';
ArtifactDataPath='D:\analysis\2008-08-27-4\data014';
ArtifactSubtraction=0;
FileName='D:\analysis\2008-08-27-4\data005\protected\clusters005'

Pattern0=317;
Movie0=73;
ps=Pattern0-97;

p1=[ps+1:ps+6];
p2=p1+12;
p3=p1+24;
p4=p1+36;
p5=p1+48;
p6=p1+60;
p7=p1+72;
p8=p1+84;

p=[p1 p1+12 p1+24 p1+36 p1+48 p1+60 p1+72 p1+84];
%eff=p;

Angles=[330 270 210 150 30 90]-30;
Angles=[270 210 330 150 90 30]-210;
for i=1:length(Angles);
    if Angles(i)<0
        Angles(i)=Angles(i)+360;
    end
end

Movies=[Movie0-2 Movie0-1 Movie0 Movie0+1];
for i=1:0
    MovieNumber=Movies(i);
    for j=1:12
        PatternNumber=p((i-1)*12+j);
        %[DataTraces,ArtifactDataTraces,Channels]=NS_ReadPreprocessedData(DataPath,ArtifactDataPath,ArtifactSubtraction,PatternNumber,MovieNumber,100);
        WaveformTypes=NS_ReadClusterFile(FileName,MovieNumber,PatternNumber);
        a=100-length(find(WaveformTypes==1));
        eff((i-1)*12+j)=a;
    end
end

eff0=51;
FigureProperties=struct('FigureNumber',5,'Subplot',[2 3 3],'TimeRange',[0 50],'AmplitudeRange',[-2 2],'FontSize',16,'Colors',['k' 'r' 'b' 'm' 'g' 'c' 'y'],'LineWidth',1,'YLabel','current [\muA]');
for i=1:4
    index=(i-1)*12
    eff1=eff(index+1:index+6);
    eff2=eff(index+7:index+12);
    FigureProperties.FigureNumber=300+i;
    h=NS_PlotEfficacyVsAngle_2008_10_13(eff0,eff1,eff2,Angles,FigureProperties);
end
