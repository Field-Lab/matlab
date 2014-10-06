a=zeros(1,1000);
b=zeros(1,1000);

%a(1,401:600)=rand(1,200)-0.5;
%b(1,401:600)=rand(1,200)-0.5;

a(1,250)=1;
b(1,240)=1;
b(1,280)=1;

c=korelacja(b,a);
d=fastcorr2(a,b,[100 100],9);

figure(1);
plot(c);
figure(2);
plot(d);

e=c(1,900:1100);

figure(3);
plot(e);

figure(4);
plot(d-e);