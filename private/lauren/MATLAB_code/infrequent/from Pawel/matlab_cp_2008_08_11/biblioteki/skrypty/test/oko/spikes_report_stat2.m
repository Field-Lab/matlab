%cd /mnt/data4/dane/2000/2000-12-12; %linux
cd I:\dane\2000\2000-12-12 %windows

nrchns=65;

s1=1000000;
dl=3000000;
t=[1:dl];

samples=[s1 (s1+dl-1)];
channels=[6 13 20 29 36 65]; % +1 bo np. kanal 1 jest w pliku jako 2
%itd.
channel=29;
%channels=31;
name='Data009';
fs=20000;

%sprawdzono, ze zmiana czestosci graqnicznej z 0.99 na np. 0.9 jest
%praktycznie bez znaczenia, dopiero przy 0.8 jest wyraznie gorzej
filtr_dl=20;
N=2;

detekcja=[400,100,100,200,100,300];

figure(3);
hold off;
clf;
figure(4);
hold off;
clf;
figure(11);
hold off;
clf;
figure(5);
hold off;
clf;

ilosc=zeros(1,6);

kolory={'bd','gd','rd','cd','md','kd'};
for i=1:6
    read_param=struct('name',name,'header',206,'nrchns',nrchns,'channels',channels(1,i),'samples',samples);
    s=readcnst(read_param);
    signal=s(1:2:dl);
    
    [y,filtr]=oversampling(signal,N,filtr_dl,0.98);
    y1=y(1,filtr_dl+1:filtr_dl+dl);
    roznica=y1-s;
    
    detect_param=struct('prog',detekcja(i),'histereza',40,'znak',-1);
    wynik=detect2(s,detect_param);
    s_wynik=size(wynik);
    blad=zeros(3,s_wynik(2));
    en=zeros(1,s_wynik(2));
    ilosc(1,i)=s_wynik(2);
    
%* * * wykresy spike'ow * * *
    figure(11);
    %hold off;
    %clf;
    subplot(2,3,i);
    for j=1:s_wynik(2);
        spike=[wynik(1,j)-10:wynik(1,j)+50];
        blad(1,j)=min(roznica(spike));
        blad(2,j)=max(roznica(spike));
	ofset=mean(s(spike));
	ofset=0;
	blad(3,j)=energia(roznica(spike),ofset);
	en(1,j)=energia(s(spike),ofset);
        plot(s(spike));
        hold on;
    end
    
    figure(2);
    subplot(2,3,i);
    plot(abs(wynik(3,:)),max(abs(blad(1,:)),abs(blad(2,:))),'bd');

%* * * figura 3 - max. bledu vs amplituda spike'a, wspolny wykres * * *
    figure(3);
    plot(abs(wynik(3,:)),max(abs(blad(1,:)),abs(blad(2,:))),kolory{i});
    axis([0 1100 0 120]);
    grid on;
    hold on;

    figure(2);
    subplot(2,3,i);
    plot(abs(wynik(3,:)),max(abs(blad(1,:)),abs(blad(2,:))),'bd');


    figure(4);
    plot(sqrt(en(1,:)),sqrt(blad(3,:)),kolory{i});
    axis([0 600 0 20]);
    grid on;
    hold on;

end
legend(['channel ' num2str(channels(1)) ': ' num2str(ilosc(1,1)) ' spikes'],['channel ' num2str(channels(2))  ': ' num2str(ilosc(1,2)) ' spikes'],['channel ' num2str(channels(3))  ': ' num2str(ilosc(1,3)) ' spikes'],['channel ' num2str(channels(4))  ': ' num2str(ilosc(1,4)) ' spikes'],['channel ' num2str(channels(5))  ': ' num2str(ilosc(1,5)) ' spikes'],['channel ' num2str(channels(6))  ': ' num2str(ilosc(1,6)) ' spikes']);
%[z,b]=resample(signal,2,1,10);
%figure(3);
%plot(z-s);
channels=[2 3 5 6 7 9 12 13 14 16 18:25 27:33 35 37 38 39 41 50 51:54 58 59 61 63 64];
l_channels=length(channels);

clear en blad wynik wynik0;
poziom=zeros(1,l_channels);
for i=1:l_channels
    read_param=struct('name',name,'header',206,'nrchns',nrchns,'channels',channels(i)+1,'samples',samples);
    s=readcnst(read_param);
    poziom(1,i)=std(s);
    i
    figura=5;
    subplot(8,5,i);
    y=blad_vs_energy(s,figura);
    %hold on;
end