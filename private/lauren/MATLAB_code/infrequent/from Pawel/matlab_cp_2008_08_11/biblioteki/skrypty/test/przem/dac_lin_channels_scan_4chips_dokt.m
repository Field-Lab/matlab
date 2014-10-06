%clear

cd H:/pliki/nauka/stymulacja/chip/testy/2006_07_12/

a=importdata('60nA_4700k_channels_scan_board1_chip_left.dat')/4700e3;
%b=importdata('4uA_75k_channels_scan_board1_chip_left.dat');
b=importdata('250uA_1190_channels_scan_board1_chip_left.dat')/1190;
%[a1,b1]=przekladki_channels_scan_dokt(a,b);

a=importdata('60nA_4700k_channels_scan_board1_chip_right.dat')/4700e3;
%b=importdata('4uA_75k_channels_scan_board1_chip_right.dat');
b=importdata('250uA_1190_channels_scan_board1_chip_right.dat')/1190;
%[a2,b2]=przekladki_channels_scan_dokt(a,b);

a=importdata('60nA_4700k_channels_scan_board3_chip_left.dat')/4700e3;
%b=importdata('4uA_75k_channels_scan_board3_chip_left.dat');
b=importdata('250uA_1190_channels_scan_board3_chip_left.dat')/1190;
%[a3,b3]=przekladki_channels_scan_dokt(a,b);

a=importdata('60nA_4700k_channels_scan_board3_chip_right.dat')/4700e3;
b=importdata('250uA_1190_channels_scan_board3_chip_right.dat')/1190;
%[a4,b4]=przekladki_channels_scan_dokt(a,b);

a1(1,21)=mean(a1(1,1:20));
a2(1,21)=mean(a2(1,1:20));
a3(1,21)=mean(a3(1,1:20));
a4(1,21)=mean(a4(1,1:20));

b1(1,21)=mean(b1(1,1:20));
b2(1,21)=mean(b2(1,1:20));
b3(1,21)=mean(b3(1,1:20));
b4(1,21)=mean(b4(1,1:20));

c1=b1./a1;
c2=b2./a2;
c3=b3./a3;
c4=b4./a4;

figure(1)
s=[1:32];
h1=plot(s,a1*1e9,s,a2*1e9,s,a3*1e9,s,a4*1e9);
set(h1,'LineWidth',2);
set(h1,'MarkerSize',8)
grid on
h1=gca;
set(h1,'LineWidth',2);
set(h1,'FontSize',18);
set(h1,'XLim',[0 33]);
set(h1,'YLim',[0.5 0.56]);
xlabel('N');
ylabel('wzmocnienie [\nA]');

figure(2)
h1=plot(s,b1*1e6,'bd-',s,b2*1e6,'rd-',s,b3*1e6,'kd-',s,b4*1e6,'cd-');
set(h1(4),'Color',[0.1 0.7 0.1])
set(h1,'LineWidth',2);
set(h1,'MarkerSize',6)
grid on
h1=gca;
set(h1,'LineWidth',2);
set(h1,'FontSize',18);
set(h1,'XLim',[0 33]);
xlabel('N');
ylabel('wzmocnienie [\muA]');

figure(3)
h1=plot(s,c1,s,c2,s,c3,s,c4);
set(h1,'LineWidth',2);
set(h1,'MarkerSize',8)
grid on
h1=gca;
set(h1,'LineWidth',2);
set(h1,'FontSize',18);
set(h1,'XLim',[0 33]);
xlabel('N');
ylabel('wzmocnienie');

a=[a1 a2 a3 a4];
b=[b1 b2 b3 b4];
c=[c1 c2 c3 c4];
std(a1)/mean(a1)*100
std(a2)/mean(a2)*100
std(a3)/mean(a3)*100
std(a4)/mean(a4)*100

std(b1)/mean(b1)*100
std(b2)/mean(b2)*100
std(b3)/mean(b3)*100
std(b4)/mean(b4)*100

std(c1)/mean(c1)*100
std(c2)/mean(c2)*100
std(c3)/mean(c3)*100
std(c4)/mean(c4)*100

std(a)/mean(a)*100
std(b)/mean(b)*100
std(c)/mean(c)*100
%plot(s,slope1,'bd-',s,slope2,'kd-',s,slope3,'rd-',s,slope4,'cd-');