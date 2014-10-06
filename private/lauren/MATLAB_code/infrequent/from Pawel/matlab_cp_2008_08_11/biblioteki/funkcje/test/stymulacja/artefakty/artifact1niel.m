function [Uc,Ir,Ic]=artifact1niel(I,t,V0,const);
%Daje przebieg napiecia na Ce, model nieliniowy, ale z pojemnoscia zamiast
%CPE. Przebieg z modelu analitycznego (nie symulacje).
%I - prad
%t - czas
%V0 - warunek poczatkowy
%const - struktura zawierajaca: I0 (typowo 9nA), Ce (1.5nF), Ut (26mV);

I0=const.I0;
Ce=const.Ce;
Ut=const.Ut;

a=I/Ce;
b=I0/Ce;
g=1/Ut;

delta=sqrt(a*a+4*b*b);

z=exp(V0/Ut);
%c1=-2/g/a*atanh((2*b*z-a)/a);
c1=-2/g/delta*atanh((2*b*z-a)/delta);

dt=t(1,2)-t(1,1);
t=t-c1;
%Uc=Ut*log(a/2/b*(1+tanh(0.5*t*g*a)));
Uc=Ut*log(1/2/b*(a+delta*tanh(0.5*t*g*delta)));

for i=2:length(t)-1
    yprim(i-1)=(Uc(1,i+1)-Uc(1,i-1))/2/dt;
end
Ic=yprim*Ce;
y1=Uc(1,2:length(Uc)-1);
Ir=I0*(exp(y1/Ut)-exp(-y1/Ut));