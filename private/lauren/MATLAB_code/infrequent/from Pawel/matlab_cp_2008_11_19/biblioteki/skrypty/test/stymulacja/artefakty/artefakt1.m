step=1e-6; %wazne!! - zmiana krok powoduje, ze wartosci artefaktu trzeba przeskalowac!
N=25000;
t=[1:N]*step; %czas w us
a1=-1;
a2=-1;
dlugosc=5000;
T=dlugosc*step; %szerokosc impulsu (jednej czesci) w us
%t=t+2*T;

art=t;
R=0.0001;

n=0.8;
t1=[1:dlugosc]*step;
art(1,1:dlugosc)=a1*t1.^n-R;
t1=[dlugosc+1:2*dlugosc]*step;
art(1,dlugosc+1:2*dlugosc)=a1*(t1.^n-(t1-T).^n)-a2*(t1-T).^n+R; %bo skonczyla sie pierwsza czesc impuslu i jednoczesnie zaczela druga!
t1=[2*dlugosc+1:N]*step;
art(1,2*dlugosc+1:N)=a1*((t1).^n-(t1-T).^n)-a2*((t1-T).^n-(t1-2*T).^n);
art=art*20;
%art1=art(1,2001:100000);
figure(1)
plot(art)
%axis([1 1000 -0.01 0.01]);
grid on
