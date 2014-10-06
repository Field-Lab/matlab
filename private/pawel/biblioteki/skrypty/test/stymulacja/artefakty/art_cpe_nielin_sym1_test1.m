clear

Y=1e-8;
beta=0.8;
Rs=20000;
I0=1e-9;
alfa=0.5;
N=1;

model_nielin=struct('Y',Y,'beta',beta,'Rs',Rs,'I0',I0,'alfa',alfa,'N',N);
[V,Icpe,Ire]=art_cpe_nielin_sym1(model_nielin,1e-6,1e-4,1.2e-3,1e-6);

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