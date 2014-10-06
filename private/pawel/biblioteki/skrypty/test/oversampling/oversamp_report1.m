filtr_dl=20;
figura1=55;
figura2=56;

channels=zeros(2,4);

cd I:\dane\2004\12maja\ForPawel\2003-08-06-0-016;
s1=importdata('electrode5.txt');
s1=s1-mean(s1);
channels(1,1)=9;
channels(2,1)=9; %do poprawki
s2=importdata('electrode399.txt');
s2=s2-mean(s2);
channels(1,2)=63;
channels(2,2)=100;
cd I:\dane\2004\12maja\ForPawel\2003-09-19-0-002;
s3=importdata('electrode203.txt');
s3=s3-mean(s3);
channels(1,3)=1; %nie! 1!
channels(2,3)=1; %do poprawki;
s4=importdata('electrode509.txt');
s4=s4-mean(s4);
channels(1,4)=77;
channels(2,4)=28;

signals=[s1 s2 s3 s4]';
dl=length(s1);
clear s1 s2 s3 s4;
ts=[1:dl]/20; %czas w milisekundach
marg_left=15;
marg_right=20;
ts=[0:marg_right+marg_left-1]/20; %czas w milisekundach

figure(figura);
for i=1:4
    s=signals(i,:);
    y=blad_vs_energy(s,28);
    a=find(y(3,:)>8); %nie zmieniac progu!!
    la=length(a)
    wskazniki=y(1,a);
    
    signal=s(1:2:dl); 
    [y0,filtr]=oversampling(signal,2,filtr_dl,0.98);
    y1=y0(1,filtr_dl+1:filtr_dl+dl);
    roznica=y1-s;
    
    sp1_start=wskazniki(channels(1,i))-marg_left;
    sp1_stop=wskazniki(channels(1,i))+marg_right-1;
    spike_wspol=[sp1_start:sp1_stop];
    spike=s(spike_wspol);
    spike=spike-mean(spike);
    subplot(3,4,i);
    plot(ts,spike,'b-',ts,y1(spike_wspol),'r-');
    axis([0 1.7 -200 200]);
    grid on;
    
    subplot(3,4,i+4);
    plot(ts,spike,'bd-',ts,y1(spike_wspol),'rd-');
    axis([0.4 0.8 -200 200]);
    grid on;
    
    subplot(3,4,i+8);
    fs=fft(spike,500);
    f=[0:499]/500*20000;
    plot(f,abs(fs));
    grid on;
    axis([0 10000 0 1000])
end