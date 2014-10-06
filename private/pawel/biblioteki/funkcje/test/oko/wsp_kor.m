function w=wsp_kor(x,y);

x1=x-mean(x);
y1=y-mean(y);

xnorma=mean(x1.^2);
ynorma=mean(y1.^2);

w=mean(x1.*y1)/sqrt(xnorma*ynorma);
