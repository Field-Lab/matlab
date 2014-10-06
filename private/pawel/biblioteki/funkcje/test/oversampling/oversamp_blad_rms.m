function s=oversamp_blad_rms(filename,detect_param,oversamp_param);
%Dane w postaci pliku tekstowego z danymi dla jendego kanalu!!!
N=oversamp_param.N;
filtr_dl=oversamp_param.filtr_dl;
filtr_gr=oversamp_param_filtr_gr;

s=importdata(filename)';
s=s-mean(s);
dl=length(s);
signal=s(1:2:dl);
    
[y,filtr]=oversampling(signal,N,filtr_dl,0.98);
y1=y(1,filtr_dl+1:filtr_dl+dl);
roznica=y1-s;
    
detect_param=struct('prog',prog,'histereza',histereza,'znak',znak);
wynik=detect2(s(1+left_marg:dl-right_marg),detect_param);
wynik=wynik+left_marg;
s_wynik=size(wynik);
   
blad=zeros(3,s_wynik(2));
en=zeros(1,s_wynik(2));
s=zeros(3,s_wynik(2));
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
s(1,:)=wynik(1,:);
s(2,:)=energia(1,:);
s(3,:)=blad(3,:);