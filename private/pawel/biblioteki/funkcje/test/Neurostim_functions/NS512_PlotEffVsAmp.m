function [h,amp]=NS512_PlotEffVsAmp(A,ClusterIndexes);

SCI=size(ClusterIndexes);

for i=1:SCI(2)
    a=find(ClusterIndexes(:,i,:)==2);
    if a>0
        pattern=i;        
    end
end

for i=1:SCI(1)
    b=ClusterIndexes(i,pattern,:);
    a=find(b==2);
    eff(i)=length(a);
end
l=find(eff>20);
amp=l(1);
%figure(1);
size(A)
size(eff)
h=plot(A,eff(1:55),'bd-');
axis([0.5 4.2 -10 110]);
grid on