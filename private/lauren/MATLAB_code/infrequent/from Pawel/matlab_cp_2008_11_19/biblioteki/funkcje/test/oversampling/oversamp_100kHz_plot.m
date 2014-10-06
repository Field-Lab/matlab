function n=oversamp_100kHz_plot(s,N,lgnd,style);
%style: 0 - rysuje cztery wykresy (raw, 10kHz, raw-N*20kHz, 10kHz-N*20kHz)
%- linie plus punkty; 1 - rysuje tylko dwa wykresy (obydwa N*20kHz), same
%linie.
%Nadprobkowanie:
fp=20;
t0=[0:length(s)-1]/fp;
czest_gr=0.95;

filtr_dl=36;
%N=10;
[y,filtr]=oversampling(s,N,filtr_dl,czest_gr);
y0=y(1,1+filtr_dl/2*N:length(y)-filtr_dl/2*N+N-1);
clear y;
t1=[0:length(y0)-1]/(N*fp);
y0a=y0(1,[length(y0)/2-20:length(y0)/2+100]);
a=find(y0a==min(y0a));
a=a+length(y0)/2-20;
n=t1(a);

signal=s(1:2:length(s));
t2=[0:length(signal)-1]/fp*2;
N=N*2;
[y1,filtr]=oversampling(signal,N,filtr_dl,czest_gr);
y2=y1(1,1+filtr_dl/2*N:length(y1)-filtr_dl/2*N+N-1);
clear y1;
t3=[0:length(y2)-1]/(N*fp)*2;

if style==0
    plot(t0,s,'bd-',t1,y0,'gd-',t2,signal,'dr-',t3,y2,'kd-');
else
    plot(t1,y0,'g-',t3,y2,'k-');
end
if lgnd==1
    if style==0
        legend('raw 20kHz','20kHz->200kHz','undersampled 10kHz','10kHz->200kHz');
    else
        legend('20kHz->200kHz','10kHz->200kHz');
    end
end
grid on;