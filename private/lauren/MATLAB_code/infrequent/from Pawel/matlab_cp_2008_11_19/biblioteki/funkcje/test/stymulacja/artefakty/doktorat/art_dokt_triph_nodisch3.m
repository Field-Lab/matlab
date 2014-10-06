function w1=art_dokt_triph_nodisch3(model,Ampl,t_disch,T,dt,fig_nr,zakres);
t0=0.1e-3;
t1=1e-4; %czas trwania pierwszej czesci impulsu
t2=1e-4; %t2 - czas trwanai drugiej czesci
t3=1e-4;

t=[dt:dt:T]*1000;

ws=[0.96 0.98 1 1.02 1.04];
for j=1:5                    
        %t_disch=200e-6;
        A1=0.6*Ampl;
        A2=-Ampl;%*ws(1,j)
        A3=-A1-A2;        
        [Vcourse,Icpe,Ire]=art_cpe_nielin_triph_no_disch(model,A1,t1,A2,t2,A3,t3,1e-4,T,dt,0);        
        w1(j,:)=Vcourse;
        size(w1);
end
figure(fig_nr);
w1=w1*1000;
zakres=zakres*1000;
a=plot(t,w1(1,:),t,w1(2,:),t,w1(3,:),t,w1(4,:),t,w1(5,:));

set(a,'LineWidth',2);
axis([0 T*1000 -zakres zakres]);
grid on;      
h=gca;
set(h,'LineWidth',2);
set(h,'FontSize',22);
set(h,'XTick',[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1]);
a=legend('I_{neg}=0.96\muA','I_{neg}=0.98\muA','I_{neg}=1.0\muA','I_{neg}=1.02\muA','I_{neg}=1.04\muA');
set(a,'LineWidth',1);

xlabel('t [ms]');
ylabel('V_{np} [mV]');