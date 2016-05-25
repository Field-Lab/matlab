clear

NS_GlobalConstants=NS_GenerateGlobalConstants(500);
cd D:\Home\Data\retina\LaurenGrosberg;
Movienumber=2;
ElectrodeNumber=11;

MovieData=NS_MovieData('005',2,NS_GlobalConstants);
break
for i=1:16239
    c1(i)=patterns_out(i).channel;
end

plot(c1,'bd-')

p1=find(c1==11)
%for j=1:length(p1)
%   a1=