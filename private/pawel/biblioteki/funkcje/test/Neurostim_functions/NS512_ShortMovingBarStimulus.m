function [electrodes,Array,Times]=NS512_LongMovingBarStimulus(DelayInMs,TimeShiftInMs,Columns);
% Be careful!! - thisis tricky, because there is 16 rows on the 512 array,
% but they are shifted by 30 microns. So each column has only 8 electrodes,
% and separation between columns is 30 microns. Here, we understand one
% column as TWO columns which equal to 16 electrodes.

ArrayClockwiseRotation=270;

TimeShift=TimeShiftInMs*20;
Delay=DelayInMs*20;

electrodes=[1:512];
Array=zeros(512,length(Columns)); %16 different patterns, each of the 512 electrodes is included in one pattern

electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(500);
x=zeros(1,512);
y=zeros(1,512);
for i=1:512
    x(i)=electrodeMap.getXPosition(i);
    y(i)=electrodeMap.getYPosition(i);
end

[a,b]=hist(x,20000);
c=find(a);
XDiff=round(b(c)); %all different x coordinates that exist in the electrode map - 64 values

[a,b]=hist(y,20000);
c=find(a);
YDiff=round(b(c));

for i=1:length(Columns)
    index=Columns(i);
    el1=find(x==XDiff(index)); %all the electrodes that have this specific value of y coordinate
    el2=find(x==XDiff(index+1));
    Array(el1,i)=1;
    Array(el2,i)=1;
    Times(i)=Delay+TimeShift*(i-1);
end