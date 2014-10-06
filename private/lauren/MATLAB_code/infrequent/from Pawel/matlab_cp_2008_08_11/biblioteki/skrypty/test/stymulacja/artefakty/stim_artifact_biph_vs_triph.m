%clear;
f=[5:1:200000];

c=36e+7;
Rs=60000;
beta=0.86;
%beta=1;
I0=25e-10;

zcpe=c*(i*2*pi*f).^(-beta);

t0=0.1e-3;
t1=1e-4; %czas trwania pierwszej czesci impulsu
t2=1e-4; %t2 - czas trwanai drugiej czesci
t3=1e-4;

Ampl=10e-7;
t_delay=-0.2e-5;
T=2.5e-3; %T - czas symulacji
dt=2e-6; %dt - krok czasowy

figure(7)
a=loglog(f,abs(zcpe)+Rs);
set(a,'LineWidth',2);
grid on;
axis([10 200000 1e4 2e7]);
h=gca;
set(h,'FontSize',24);
xlabel('f [Hz]');
ylabel('impedance module [k\Omega]')

%break;

%cd H:\pliki\nauka\doktorat\obrazki\rozdzial4;
% 1. Impuls dwufazowy, model liniowy, bez rozladowania
model=struct('Y',1/c,'beta',beta,'Rs',Rs,'I0',I0/1e6,'alfa',0.5,'N',1); 
figura=221;
figure(figura);
%w1=art_dokt_biph_nodisch3(model,Ampl,200e-6,T,dt,figura,2e-3);

% 2. Impuls dwufazowy, model liniowy, z rozl.
model=struct('Y',1/c,'beta',beta,'Rs',Rs,'I0',I0,'alfa',0.5,'N',1); 
figura=312;
figure(figura);
%w1=art_dokt_biph_disch(model,Ampl,200e-6,T,dt,figura,2e-3);
%print(figura,'-dpng','ArtBiphLinDisch80kOhm.png');

% 2. Impuls dwufazowy, model liniowy, z rozl.
model=struct('Y',1/c,'beta',beta,'Rs',Rs,'I0',5e-15,'alfa',0.5,'N',1); 
figura=22;
figure(figura);
%w1=art_dokt_biph_disch(model,10e-7,200e-6,T,dt,figura,150e-3);
%print(figura,'-dpng','ArtBiphLinDisch80kOhm.png');

% 3. Impuls trzyfazowy, model liniowy, bez rozl.
figura=3;
figure(figura);
model=struct('Y',1/c,'beta',beta,'Rs',Rs,'I0',5e-15,'alfa',0.5,'N',1); 
w1=art_dokt_triph_nodisch(model,1e-6,200e-6,T,dt,figura,2e-3);%
%print(figura,'-dpng','ArtTriphLinNodisch.png');

% 4. Jak wyzej, beta=0.82;
model=struct('Y',1/c,'beta',0.82,'Rs',Rs,'I0',5e-15,'alfa',0.5,'N',1); 
figura=4;
%figure(figura);
%w1=art_dokt_triph_nodisch(model,1e-6,200e-6,T,dt,figura,4e-3);
%print(figura,'-dpng','ArtTriphLinNodischBeta082.png');

% 5. Jak wyzej, beta=0.90;
model=struct('Y',1/c,'beta',0.9,'Rs',Rs,'I0',5e-15,'alfa',0.5,'N',1); 
figura=5;
%figure(figura);
%w1=art_dokt_triph_nodisch(model,1e-6,200e-6,T,dt,5,4e-3);
%print(figura,'-dpng','ArtTriphLinNodischBeta090.png');

% 6. Impuls dwufazowy, model nieliniowy, bez rozladowania
%figura=6;
%figure(figura);
model=struct('Y',1/c,'beta',beta,'Rs',Rs,'I0',25e-9,'alfa',0.5,'N',1); 
%w1=art_dokt_biph_nodisch2(model,0.8e-6,200e-6,T,dt,figura,20e-3);
%print(figura,'-dpng','ArtBiphNielinDisch80kOhm1uA.png');

% 7. Impuls dwufazowy, model nieliniowy, z rozl.
figura=27;
figure(figura);
%model=struct('Y',1/c,'beta',beta,'Rs',80000,'I0',25e-9,'alfa',0.5,'N',1); 
%w1=art_dokt_biph_disch(model,4e-6,200e-6,T,dt,figura,40e-3);
%print(figura,'-dpng','ArtBiphNielinDisch80kOhm4uA.png');

% 8. Impuls trzyfazowy, model nieliniowy, bez rozl.
figura=28;
figure(figura);
%model=struct('Y',1/c,'beta',beta,'Rs',Rs,'I0',25e-9,'alfa',0.5,'N',1); 
model=struct('Y',1/c,'beta',beta,'Rs',Rs,'I0',I0/2,'alfa',0.5,'N',2); 
%w1=art_dokt_triph_nodisch(model,1e-6,200e-6,T,dt,figura,4e-3);
%print(8,'-dpng','ArtTriphNolinNodisch5uAWar.png');

% 8. Impuls trzyfazowy, model nieliniowy, z rozl.
figura=9;
%figure(figura);
%model=struct('Y',1/c,'beta',beta,'Rs',80000,'I0',25e-9,'alfa',0.5,'N',1); 
%w1=art_dokt_triph_disch3(model,5e-6,200e-6,T,dt,figura,4e-3);
%print(9,'-dpng','ArtTriphNolinDisch80kOhm5uAWar.png');

figura=19;
%figure(figura);
model=struct('Y',1/c,'beta',beta,'Rs',80000,'I0',25e-9,'alfa',0.5,'N',1); 
%w1=art_dokt_triph_disch(model,5e-6,200e-6,T,dt,figura,4e-3);
%print(19,'-dpng','ArtTriphNolinDisch400kOhm5uA.png');