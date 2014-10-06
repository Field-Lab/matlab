function [y,psdf,sigma]=crosstalk1(filename,fs,f0,channel0,wykladnik,length);
%Funkcja zwraca wektor wspol. fourierowskich dla czestotliwosci
%f0 (po wszystkich kanalach) oraz komplet funkcji spdf dla wszystkich
%kanalow, a takze wartosc sredniokwadratowa (sigme) z calego sygnalu (dla
%wszystkich kanalow).
%Rozdzielczosc czest. obu danych wynikowych jest identyczna i zdefiniowana
%przez wykladni oraz fs: df=fs/(2^wykladnik). 
%Dlugosc danych wynosi typowo 600000, 
%dlugosc fragmentu danych poddawanego fft  typowo 16384. Mozna wiec wykonac
%kilkadziesiat usrednien.
%Wektor y liczy sie nastepujaco: dla jednej iteracji wyznacza sie tranf.
%fourier. dla wszystkich kanalow, nastepnie dla czest. f0 (dokladnie - dla
%najblizszej czest sposrod bazowych czest fft) dokonuje sie korekty fazy.
%Faza sygnalu w kanale channel0 (gdzie podawany byl sygnal, zatem poziom
%jest wysoki i faza nie jest zaburzona przez szumy) jest odejmowana od fazy
%wspol fourier. dla tej czest. we wszystkich kanalach. Tak utworzon wetkor
%y jest zapamietywany. W nastepnej iteracji znow obliczane jest fft,
%forygowana faza, i wynik dodawany (bezp. wspol fourier., nie wartosci
%bezwzgledne) do poprzedniego. Dzieki korekcie fazy jest to "usrednianie
%synchroniczne). Wartosci sa w jednostkach LSB.
%Wartosci spdf liczone sa na bazie tych samuch transformat fft co wektor y.
%Usredniane sa jednak oczywiscie moduly tranformat a nie wspol. zespolone
%fft. Wartosci sa wyrazone w LSB^2/Hz.
%length=600000;
mnoznik=1.5;
header=108;

N=2^wykladnik;
steps=floor(length/N)
%steps=3;
f_nr=round(N*f0/fs)+1 % nr czestotliwosci f0 w tr. fouriera
df=fs/N; %rozdzielczosc analizy fft, wazne dla obliczenia SPDF

fid = fopen(filename,'r');
fseek(fid,header+2,-1);

F=zeros(512*mnoznik,N);
fsig=zeros(512,N);

y=zeros(1,512);
dane=zeros(1,512);
sigma=zeros(1,512);
psdf=zeros(512,N);

for i=1:steps
    i
    %1. czytanie danych - N probek
    for j=1:N
        F(:,j)=fread(fid,512*mnoznik,'ubit8',0);
        fseek(fid,2,0);
    end
    
    %2. dekodowanie kanalow (N probek) i fft
    for j=1:256    
        offset=2*mnoznik*(j-1);
        b1=F(offset+1,:);
        b2=F(offset+2,:);
        b3=F(offset+3,:);

        s1=b1*16+floor(b2/16)-2048;
        s2=(b2-floor(b2/16)*16)*256+b3-2048;
        
        s1=s1-mean(s1);
        s2=s2-mean(s2);
    
        %fsig(2*j-1,:)=fft(s1);
        %fsig(2*j,:)=fft(s2);
        
        f1=fft(s1)/N*2;
        f2=fft(s2)/N*2;
        
        % 2a) wektor wspol. dla wyroznionej czest. (f0) po wszystkich
        % kanalach
        dane(1,2*j-1)=f1(1,f_nr);
        dane(1,2*j)=f2(1,f_nr);
        
        sigma(1,2*j-1)=sigma(1,2*j-1)+std(s1)^2;
        sigma(1,2*j)=sigma(1,2*j)+std(s2)^2;
        
        % 2b) gestosc widmowa:
        % norm. transf. fft
        %f1=f1/N*2;
        %kwadrat modulu
        f1=f1.*conj(f1);
        %uwzglednienie kwantu czest. dla fft
        f1=abs(f1)/df;
        %wartosc skuteczna sinusa to A/sqrt(2), wiec:
        f1=f1/2;               
        psdf(2*j-1,:)=psdf(2*j-1,:)+f1;
        
        %f2=f2/N*2;
        f2=f2.*f2;
        f2=abs(f2)/df;
        f2=f2/2;        
        psdf(2*j,:)=psdf(2*j,:)+f2;
    end
       
    %3. korekta fazy        
    wspol=dane(1,channel0);
    korekta=conj(wspol)/abs(wspol);
    dane=dane*korekta;
    y=y+dane;
end
y=y/steps;
psdf=psdf/steps;
sigma=sqrt(sigma/steps)