clear;
cm=1; %uF na cm kw.
cm=cm/100; %w jedn. SI
rho_i=70; %ohm*cm
rho_i=rho_i/100; %w jedn. SI
rho_e=200; %ohm*cm

d=5e-6; %srednica aksonu;
dx=1e-6; %kwantowanie dlugosci aksonu
cm=cm*pi*dx*d

rho_i=rho_i*4*dx/(pi*d*d)

el_pos=500;
el_odl=30e-6;
N=1000;
Vi=zeros(1,N); %intracellular potential along the fiber;
Ve=zeros(1,N); %extracellular potential along the fiber;
Ve(1,el_pos)=0.1;

for i=1:N
    sq=sqrt((el_pos-i)^2*dx^2+el_odl^2);
    Ve(i)=1e-5/sq;
end
figure(2)
plot(Ve)

figure(1)
dt=1e-5;
T=1e-3;

t=[dt:dt:T];

Vi=Ve; % stan poczatkowy - po wlaczeniu pola zewnetrznego
Vi_next=Vi;
for i=1:50 %length(t)
    for j=2:N-1
        dV=(2*Vi(j)-Vi(j-1)-Vi(j+1))/(rho_i*cm);
        Vi_next(j)=Vi(j)-dV*dt;
    end
    i
    plot(Vi(1,200:800))
    pause(1)
    Vi=Vi_next;
end