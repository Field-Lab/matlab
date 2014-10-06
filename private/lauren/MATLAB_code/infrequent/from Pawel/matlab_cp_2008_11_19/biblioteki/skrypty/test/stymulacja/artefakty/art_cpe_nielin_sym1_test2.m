clear

Y=1e-8;
beta=0.8;
Rs=20000;
I0=1e-10;
alfa=0.5;
N=2;

model_nielin=struct('Y',Y,'beta',beta,'Rs',Rs,'I0',I0,'alfa',alfa,'N',N);

A=2e-6;
c1=0.5;
A1=c1*A
A2=-A
A3=-A1-A2

%[V1,Icpe1,Ire1]=art_cpe_nielin_sym1(model_nielin,A1,1e-4,1.2e-3,1e-6);
%[V2,Icpe2,Ire2]=art_cpe_nielin_sym1(model_nielin,A2,1e-4,1.1e-3,1e-6);
%[V3,Icpe3,Ire3]=art_cpe_nielin_sym1(model_nielin,A3,1e-4,1e-3,1e-6);

%V2=[zeros(1,length(V1)-length(V2)) V2];
%V3=[zeros(1,length(V1)-length(V3)) V3];

t1=1e-4;
t2=t1;
t3=t1;
T=1e-3;
dt=1e-6;
[V,Icpe,Ire]=art_cpe_nielin_sym1_triphas(model_nielin,A1,t1,A2,t2,A3,t3,T,dt);

figure(2)
subplot(2,2,1);
plot(V)
grid on;
subplot(2,2,2);
plot(Icpe)
grid on;
subplot(2,2,3);
plot(Ire)
grid on;
subplot(2,2,4);
plot(V,Ire)
grid on;