function [rhist rbin]=histc2(x,y,bins)
%bins y based on x
%vector only

rhist=zeros(length(bins)-1,1);
rbin=zeros(length(x),1);

for i=1:length(bins)-1
    t=find(x>=bins(i) & x<bins(i+1));
    rhist(i)=sum(y(t));
    rbin(t)=i;
end



