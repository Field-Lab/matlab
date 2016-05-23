function [c,smax,sigma]=ConvPeakMax(s,w);

c=conv(s,w,'valid');
%c=c0(length(w):length(s))/length(w);

%figure(12)
%plot(c)

cmax=max(c);
cmax_ind=mean(find(c==max(c)));

i1=find(sign(c-cmax/2)>0);
i2=find(diff(i1)>1);

if i2
    i2prim=sort([min(i1) i1(i2) i1(i2+1) max(i1)])
    for i3=1:length(i2prim)-1
        i3;
        i2prim(i3);
        i2prim(i3+1);
        if i2prim(i3)<cmax_ind & cmax_ind<i2prim(i3+1)
            %i1(i2prim(i3))
            %i1(i2prim(i3+1))
            center=(i2prim(i3)+i2prim(i3+1))/2;
            sigma=(i2prim(i3+1)-i2prim(i3))/2;
        end
        if i2prim(i3)<cmax_ind & cmax_ind==i2prim(i3+1)
            center=i2prim(i3+1);
            sigma=(i2prim(i3+1)-i2prim(i3))/2;
        end
        if i2prim(i3)==cmax_ind & cmax_ind<i2prim(i3+1)
            center=i2prim(i3);
            sigma=(i2prim(i3+1)-i2prim(i3))/2;
        end
    end
else
    center=mean(i1);
    sigma=(max(i1)-min(i1))/2;
end
smax=center+(length(w)+1)/2-1;
%smax=mean(cmax+(length(w)+1)/2-1);