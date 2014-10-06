cd /home/pawel/pliki/nauka/neuroplat;

ilf=20;
ihf=13;
igain=24;

a=importdata('pasmo_ch15_DAC_cal1_ilf20_ihf13_gain24_restfloating.dat');
bm=a(17:21,:);
bm1=bm(:,2:29)*1e6/135;
b=bm(1,:);
c=importdata('czest.dat');
c1=c(1,1:length(c));
b1=b(1,2:29)*1e6/135;

% 1. Pomiar z zewnetrznym generatorem:
% a) czestotliwosci - wg kolejnosci pomiaru:
f=[11 33 55 19 57 95 38 114 190 76 228 380 156 468 780 306 918 1530 612 1836 3060 1220 3660 6100 2440 7320];

% b) wspolczynniki wg pomiaru sinusem:
g1=[24 89 139 49 143 195 102 211 242 173 249 259 232 260 256 257 252 231 259 219 173 242 155 105 195 89]/1000;

% c) wspolczynniki wg pomiaru prostokatem:
g2=[23 87 136 48 140 191 100 206 233 169 239 246 229 252 247 254 249 229 258 218 173 241 154 104 194 88]/1000;

dane0=[f' g1' g2']';
dane=sort(dane0,2);

[f0,i]=sort(f);
g10=g1(i)*1e6/234;
g20=g2(i)*1e6/234;

mnoznik=234/294; %z pomiaru:wysokosci piku oraz rzeczywistej amplitudy
g10=g10*mnoznik;
g20=g20*mnoznik;
b1=b1*mnoznik;

figure(2);
loglog(f0,g10,'bd-',f0,g20,'r*-',c1,b1,'k+-');

%loglog(c1,b1,'bd-',f,g,'r*-');
grid on;
fontsize=18;
axis([10 10000 50 1000]);
h1=gca;
set(h1,'FontSize',18);
legend('external sinus','external square','internal calibration');
xlabel('frequency');
set(h1,'FontSize',18);
ylabel('gain');
set(h1,'FontSize',18);


figure(3);
bm1=bm(:,2:29)*1e6/135*mnoznik;
loglog(c1,bm1(1,:),c1,bm1(2,:),c1,bm1(3,:),c1,bm1(4,:),c1,bm1(5,:));

%loglog(c1,b1,'bd-',f,g,'r*-');
grid on;
fontsize=18;
axis([10 10000 50 1000]);
h1=gca;
set(h1,'FontSize',18);
%legend('external sinus','external square','internal calibration');
h1=xlabel('frequency');
set(h1,'FontSize',18);
h1=ylabel('gain');
set(h1,'FontSize',18);


figure(3)
clf;
subplot(2,2,1);

max11=max(g10);
czest1=[20:0.5:5000];

maxg10=max(g10)
maxg20=max(g20)
maxb1=max(b1)

s10=spline(f0,g10,czest1);
s20=spline(f0,g20,czest1);
s1=spline(c1,b1,czest1);

a=find(s10>maxg10/sqrt(2));
flow=czest1(min(a))
fhigh=czest1(max(a))

a=find(s20>maxg20/sqrt(2));
flow=czest1(min(a))
fhigh=czest1(max(a))

a=find(s1>maxb1/sqrt(2));
flow=czest1(min(a))
fhigh=czest1(max(a))

