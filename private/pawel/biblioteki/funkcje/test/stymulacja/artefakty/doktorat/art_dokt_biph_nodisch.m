function w1=art_dokt_biph_nodisch(model,Ampl,t_disch,T,dt,fig_nr,zakres);
t0=0.1e-3;
t1=1e-4; %czas trwania pierwszej czesci impulsu
t2=1e-4; %t2 - czas trwanai drugiej czesci
t3=1e-4;

t=[dt:dt:T]*1000;

for j=1:1                    
        t_disch=200e-6;
        A1=0;%(0.4+0.05*j)*Ampl
        A2=-Ampl;
        A3=-A1-A2;        
        [Vcourse,Icpe,Ire]=art_cpe_nielin_triph_no_disch(model,A1,t1,A2,t2,A3,t3,1e-4,T,dt,0);        
        w1(j,:)=Vcourse;
        size(w1);
end
figure(fig_nr);
%a=plot(t,w1(1,:),t,w1(2,:),t,w1(3,:),t,w1(4,:),t,w1(5,:),t,w1(6,:));
w1=w1*1000;
zakres=zakres*1000;

subplot(3,1,2);
v=w1(1,:)+Ampl*130000;
plot(v);

subplot(3,1,3);
a=plot(t,w1(1,:));
set(a,'LineWidth',2);
axis([0 T*1000 -zakres zakres]);
grid on;      
h=gca;
set(h,'LineWidth',2);
set(h,'FontSize',22);
set(h,'XTick',[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1]);
xlabel('t [ms]');
ylabel('V_np [mV]');