function [XDiff,YDiff]=NS512_512ElectrodeMapAllCoordinates;

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