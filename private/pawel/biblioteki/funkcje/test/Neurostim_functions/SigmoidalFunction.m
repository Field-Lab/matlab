function y=SigmoidalFunction(Xscaling,Xshift,Yscaling,Yshift,x)

y0=(1+exp((x-Xshift)/Xscaling)).^-1;
y=Yscaling*y0+Yshift;