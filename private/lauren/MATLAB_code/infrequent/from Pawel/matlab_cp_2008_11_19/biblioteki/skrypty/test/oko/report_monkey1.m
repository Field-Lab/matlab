cd I:\dane\2004\12maja\ForPawel\2003-08-06-0-016;
cd I:\dane\2004\12maja\ForPawel\2003-09-19-0-002;
cd I:\dane\2004\26maja\2003-12-08-0-001;
%cd I:\dane\2004\26maja\2004-01-13-0-001;

N=2;
filtr_dl=20;

prog=50;
histereza=40;
znak=-1;

left_marg=10;
right_marg=50;

pliki={'electrode5.txt' 'electrode10.txt' 'electrode78.txt' 'electrode167.txt' 'electrode168.txt' 'electrode190.txt' 'electrode200.txt' 'electrode353.txt' 'electrode399.txt' 'electrode506.txt'};
pliki={'electrode119.txt' 'electrode183.txt' 'electrode193.txt' 'electrode203.txt' 'electrode208.txt' 'electrode237.txt' 'electrode239.txt' 'electrode312.txt' 'electrode340.txt' 'electrode509.txt'};
pliki={'electrode45.txt' 'electrode137.txt' 'electrode148.txt' 'electrode155.txt' 'electrode169.txt' 'electrode221.txt' 'electrode298.txt' 'electrode362.txt' 'electrode416.txt' 'electrode499.txt'};
%pliki={'electrode35.txt' 'electrode130.txt' 'electrode152.txt' 'electrode313.txt' 'electrode337.txt' 'electrode370.txt' 'electrode424.txt' 'electrode471.txt' 'electrode474.txt' 'electrode492.txt'};


ilosc=zeros(1,6);

kolory={'bd','gd','rd','cd','md','kd','b+','g+','r+','c+'};
figura1=23;
figura2=24;
figure(figura1);
hold off;
clf;
figure(figura2);
hold off;
clf;
for i=1:10
    s=importdata(pliki{i})';
    s=s-mean(s);
    figure(1);
    plot(s);
    dl=length(s);
    signal=s(1:2:dl);
    
    [y,filtr]=oversampling(signal,N,filtr_dl,0.98);
    y1=y(1,filtr_dl+1:filtr_dl+dl);
    roznica=y1-s;
    
    detect_param=struct('prog',prog,'histereza',histereza,'znak',znak);
    wynik=detect2(s(1+left_marg:dl-right_marg),detect_param);
    wynik=wynik+left_marg;
    s_wynik=size(wynik)
    
    blad=zeros(3,s_wynik(2));
    en=zeros(1,s_wynik(2));
    ilosc(1,i)=s_wynik(2);

    for j=1:s_wynik(2);
        spike=[wynik(1,j)-10:wynik(1,j)+50];
        blad(1,j)=min(roznica(spike));
        blad(2,j)=max(roznica(spike));
	    ofset=mean(s(spike));
	    ofset=0;
	    blad(3,j)=energia(roznica(spike),ofset);
	    en(1,j)=energia(s(spike),ofset);
        %plot(s(spike));
        %hold on;
    end
%* * * figura 3 - max. bledu vs amplituda spike'a, wspolny wykres * * *
    figure(23);
    subplot(2,5,i);
    plot(abs(wynik(3,:)),max(abs(blad(1,:)),abs(blad(2,:))),kolory{i});
    axis([0 1100 0 120]);
    grid on;
    %hold on;

    figure(24);
    subplot(2,5,i);
    plot(sqrt(en(1,:)),sqrt(blad(3,:)),kolory{i});
    axis([0 600 0 20]);
    grid on;
    %hold on;    
end
