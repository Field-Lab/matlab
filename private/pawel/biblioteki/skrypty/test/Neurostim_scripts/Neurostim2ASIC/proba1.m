VgVth=0.1414;
Vgsigma=0.001;

c1=1e-5;
W=40;
L=0.4;

I1=c1*W/L*VgVth^2
I2=c1*W/L*(VgVth+Vgsigma)^2

Error=I2/I1