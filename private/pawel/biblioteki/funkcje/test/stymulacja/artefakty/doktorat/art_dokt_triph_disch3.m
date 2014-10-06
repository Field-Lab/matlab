function w1=art_dokt_triph_disch3(model,Ampl,t_disch,T,dt,fig_nr,zakres);
t0=0.1e-3;
t1=1e-4; %czas trwania pierwszej czesci impulsu
t2=1e-4; %t2 - czas trwanai drugiej czesci
t3=1e-4;

t=[dt:dt:T]*1000;

figure(1);
ws=[0.96 0.98 1 1.02 1.04];
for j=1:5                    
        %t_disch=200e-6;
        A1=0.65*Ampl;
        A2=-Ampl*ws(1,j)
        A3=-A1-A2;        
        [Vcourse,Icpe,Ire,Irs]=art_cpe_nielin_triph_disch(model,A1,t1,A2,t2,A3,t3,1e-4,t_disch,T,dt,0);        
        w1(j,:)=Vcourse;
        size(w1);
        %subplot(2,3,j);
        plot(Vcourse,Ire);
        hold on;
        grid on;
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
xlabel('t [ms]');
ylabel('V_{np} [mV]');
a=legend('I_{neg}=4.8\muA','I_{neg}=4.9\muA','I_{neg}=5.0\muA','I_{neg}=5.1\muA','I_{neg}=5.2\muA');
