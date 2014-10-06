L=1;
dx=0.001;
zasieg=6;
X=[-zasieg:dx:zasieg];

A0=sqrt(L^2+X.^2);
A1=sqrt(L^2+(X+dx).^2);
A2=sqrt(L^2+(X-dx).^2);

fa=-(2./A0-1./A1-1./A2)*(1/dx^2);
a=plot(X,fa)
axis([-zasieg zasieg -2 0.4])
set(a,'LineWidth',2)
set(a,'Color',[51/256,153/256,102/256])
%grid on
%h=gca;
%set(h,'LineWidth',2)
min(fa)
mean(fa)