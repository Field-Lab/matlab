clear;
c=36e+6;
f=[20:1:20000];

Rs=130000;
Re=2e9;

beta=0.86;
zcpe=c*(i*2*pi*f).^(-beta);
a=zcpe+Re;

model=struct('Ce',c,'beta',beta,'Re',Re,'Rs',20000);
t0=0.1e-3;
t1=1e-4; %czas trwania pierwszej czesci impulsu
A1=0.6e-6; %amplituda pierwszej czesci
%A1=0;
t2=1e-4; %t2 - czas trwanai drugiej czesci
t2=t1;
A2=-1e-6; %A2 - amplituda drugiej czesci
t3=1e-4;
t3=t1;
A3=-A1-A2;
%A3=A3*0.95;
t_delay=1e-5;
t_disch=2e-4;
%t_disch=0;
%A2=0;
T=1e-3; %T - czas symulacji
dt=0.4e-6; %dt - krok czasowy
t=[dt:dt:T]*1000;

for i=1:6
    beta=0.6+i*0.05
    model=struct('Ce',c,'beta',beta,'Re',Re,'Rs',5000);
    for j=1:6
        A1=(0.45+j*0.05)*1e-6;
        A3=-A1-A2;
        w1(j,:)=artifact2(model,0,t0,t1,A1,t2,A2,t3,A3,t_delay,0,T,dt);
        %w2(j,:)=artifact2(model,0,t0,t1,A1,t2,A2,t3,A3,t_delay,t_disch,T,dt);
    end
    w1=w1/(max(w1(1,:))-min(w1(1,:)))*1000;
    %w2=w2/(max(w2(1,:))-min(w2(1,:)))*2500;
    figure(11);
    subplot(3,2,i);    
    plot(t,w1(1,:),t,w1(2,:),t,w1(3,:),t,w1(4,:),t,w1(5,:),t,w1(6,:));
    %axis([0 T*1000 -50 50]);
    grid on;  
    
    %figure(3);
    %subplot(3,2,i);    
    %subplot(3,2,1);    
    %plot(t,w2(1,:),t,w2(2,:),t,w2(3,:),t,w2(4,:),t,w2(5,:),t,w2(6,:));
    %axis([0 T*1000 -5 50]);
    %grid on;  
end

T=3e-3; %T - czas symulacji
dt=0.4e-6; %dt - krok czasowy
t=[dt:dt:T]*1000;
clear w1;
clear w2;
t_disch=0.25e-3;
for i=1:6
    A1=0;
    beta=0.6+i*0.05
    model=struct('Ce',c,'beta',beta,'Re',Re,'Rs',250000);
    for j=1:6
        A3=(0.75+j*0.05)*1e-6;        
        w1(j,:)=artifact2(model,0,t0,t1,A1,t2,A2,t3,A3,t_delay,0,T,dt);
        w2(j,:)=artifact2(model,0,t0,t1,A1,t2,A2,t3,A3,t_delay,t_disch,T,dt);
    end
    w1=w1/(max(w1(1,:))-min(w1(1,:)))*1000;
    w2=w2/(max(w2(1,:))-min(w2(1,:)))*1000;
    figure(12);
    subplot(3,2,i);    
    plot(t,w1(1,:),t,w1(2,:),t,w1(3,:),t,w1(4,:),t,w1(5,:),t,w1(6,:));
    %axis([0 T*1000 -50 50]);
    grid on;  
    
    figure(14);
    subplot(3,2,i);    
    plot(t,w2(1,:),t,w2(2,:),t,w2(3,:),t,w2(4,:),t,w2(5,:),t,w2(6,:));
    %axis([0 T*1000 -5 50]);
    grid on;  
end