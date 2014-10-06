function [Vcourse,Icpe,Ire]=art_cpe_nielin_sym3(model_nielin,A1,T1,A2,T2,A3,T3,T,dt,Vpocz);
%Daje przebieg napiecia na elektrodzie w wyniku przejscie impulsu pradowego
%o amplitudzie A i czasie trwania t0. Zwraca napiecie dla czasu od poczatku
%impulsu do T. 
%model_nielin - model elektrody
%t0 - czas trwania impulsu
%T - symulowany czas
%dt - krok czasowy
F=6.02e4*1.602; %stala Avogadro razy ladunek elementarny, zredukowane wykladniki przy potegach 10
R=8.314472; %stala gazowa;
Temp=300; %temperatura

f=F/(R*Temp);

Y=model_nielin.Y;
beta=model_nielin.beta;
Rs=model_nielin.Rs;
I0=model_nielin.I0;
alfa=model_nielin.alfa; %transfer coefficient, typowo 0.5
N=model_nielin.N; %ilosc elektronow na oxydacje/redukcje

t1=[0:dt:T1-dt];
t2=[0:dt:T1+T2-dt];
t3=[0:dt:T1+T2+T3-dt];
t4=[0:dt:T-dt];
Ire=zeros(1,length(t2));
Icpe=Ire;
V=Ire;

Vpocz=0;
Vi=Vpocz;
for i=1:length(t1)    
    Ire(i)=I0*(-exp(-alfa*N*f*Vi)+exp((1-alfa)*N*f*Vi));
    Icpe(i)=A1-Ire(i);
    Vi=0;
    for j=1:i      
        Vi=Vi+Icpe(j)/(Y*gamma(beta+1))*((dt*(i-j+1))^beta-(dt*(i-j))^beta);
    end
    V(i)=Vi;
end

start=i;
for i=start+1:length(t2)    
    Ire(i)=I0*(-exp(-alfa*N*f*Vi)+exp((1-alfa)*N*f*Vi));
    Icpe(i)=A2-Ire(i);
    Vi=0;
    for j=1:i      
        Vi=Vi+Icpe(j)/(Y*gamma(beta+1))*((dt*(i-j+1))^beta-(dt*(i-j))^beta);
    end
    V(i)=Vi;
end

start=i;
for i=start+1:length(t3)    
    Ire(i)=I0*(-exp(-alfa*N*f*Vi)+exp((1-alfa)*N*f*Vi));
    Icpe(i)=A3-Ire(i);
    Vi=0;
    for j=1:i      
        Vi=Vi+Icpe(j)/(Y*gamma(beta+1))*((dt*(i-j+1))^beta-(dt*(i-j))^beta);
    end
    V(i)=Vi;
end

start=i;
for i=start+1:length(t4)    
    Ire(i)=I0*(-exp(-alfa*N*f*Vi)+exp((1-alfa)*N*f*Vi));
    Icpe(i)=-Ire(i);
    Vi=0;
    for j=1:i      
        Vi=Vi+Icpe(j)/(Y*gamma(beta+1))*((dt*(i-j+1))^beta-(dt*(i-j))^beta);
    end
    V(i)=Vi;
end

Vcourse=V;
