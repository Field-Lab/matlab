function w1=art_dokt_triph_disch(model,Ampl,t_disch,T,dt,fig_nr,zakres);
t0=0.1e-3;
t1=1e-4; %czas trwania pierwszej czesci impulsu
t2=1e-4; %t2 - czas trwanai drugiej czesci
t3=1e-4;

t=[dt:dt:T]*1000;

figure(1);
for j=1:7                    
        %t_disch=200e-6;
        A1=(0.4+0.05*j)*Ampl
        A2=-Ampl;%*0.98;
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
a=plot(t,w1(1,:),t,w1(2,:),t,w1(3,:),t,w1(4,:),t,w1(5,:),t,w1(6,:),t,w1(7,:));

set(a,'LineWidth',2);
axis([0 T*1000 -zakres zakres]);
grid on;      
h=gca;
set(h,'LineWidth',2);
set(h,'FontSize',20);
set(h,'XTick',[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1]);
xlabel('t [ms]');
ylabel('V_{np} [mV]');
a=legend('A1=0.45','A1=0.50','A1=0.55','A1=0.60','A1=0.65','A1=0.70','A1=0.75');
set(a,'LineWidth',2);