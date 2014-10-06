function w2=art_dokt_biph_disch(model,Ampl,t_disch,T,dt,fig_nr,zakres);
t0=0.1e-3;
t1=1e-4; %czas trwania pierwszej czesci impulsu
t2=1e-4; %t2 - czas trwanai drugiej czesci
t3=1e-4;

t=[dt:dt:T]*1000;

for j=1:5                    
        t_disch=300e-6*(j-1)+100e-6
        A1=0;
        A2=-Ampl;
        A3=-A1-A2;   
        %A3=A3*0.95; !!!!!
        [Vcourse,Icpe,Ire,Irs]=art_cpe_nielin_triph_disch(model,A1,t1,A2,t2,A3,t3,1e-4,t_disch,T,dt,0);        
        w2(j,:)=Irs;
        w1(j,:)=Vcourse;
        size(w1);
end
figure(fig_nr);
w1=w1*1000;
zakres=zakres*1000;
t=t-0.4;
a=plot(t,w1(1,:),t,w1(2,:),t,w1(3,:),t,w1(4,:),t,w1(5,:));
%a=plot(t,w1(1,:),t,w1(2,:),t,w1(3,:),t,w1(4,:),t,w1(5,:),t,w1(6,:));

set(a,'LineWidth',2);
axis([0 T*1000 -zakres zakres]);
grid on;      
h=gca;
set(h,'LineWidth',2);
set(h,'FontSize',20);
set(h,'XTick',[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1.1 1.2]);
a=legend('t_d=100\mus','t_d=200\mus','t_d=300\mus','t_d=400\mus','t_d=500\mus');
%a=legend('t_r=50\mus','t_r=100\mus','t_r=150\mus','t_r=200\mus','t_r=250\mus','t_r=300\mus');
set(a,'LineWidth',2);
%legend('t_r=25\mus','t_r=50\mus','t_r=75\mus','t_r=100\mus','t_r=125\mus','t_r=150\mus');
xlabel('t [ms]');
ylabel('V_{el} [mV]');
%cd H:\pliki\nauka\doktorat\obrazki\rozdzial4;
