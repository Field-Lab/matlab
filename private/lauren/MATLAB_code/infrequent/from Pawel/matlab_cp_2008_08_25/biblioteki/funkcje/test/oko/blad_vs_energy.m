function y=blad_vs_energy(s,figura);
%Sprawdzenie funkcji do oversamplingu, s - sygnal dla jednego kanalu.
%Funkcja nanosi na rysunek o numerze "figura" punkty odpowiadajace kolejnym
%spike'om, o wspolrzednych (energia_spike'a,energia_bledu). Parametry
%detekcji definiowane sa wewnatrz funckji (na poczatku).
%Do liczenia energii uzywana jest funkcja energia. 
%Funkcja zwraca rozwniez tablice 3xN, w kolejnych wierszach kolejno -
%poczatek spika (nr probki), energia spika (rms), energia bledu (rms).
%1. Definicja stalych:
prog=100;
histereza=40;
znak=-1;

N=2;
filtr_dl=16;
f_gr=0.98;

%2. Undersampling:
dl=length(s);
signal=s(1:2:dl); 

%3. Oversampling:
[y,filtr]=oversampling(signal,N,filtr_dl,f_gr);
y1=y(1,filtr_dl+1:filtr_dl+dl);
roznica=y1-s;
        
%4. Detekcja spike'ow:
detect_param=struct('prog',prog,'histereza',histereza,'znak',znak);
wynik0=detect2(s,detect_param); %"wynik" zawiera tylko numery probek, a nie probki!

%5. Ekstrakcja spike'ow: od punktu, w ktorym wykryto spike'a, bierzemy pewna
%ilosc probek "w tyl" i pewna "w przod", otrzymujac pelny ksztalt
%pojedynczego spike'a
left=10;
right=20;

wynik=wynik0; %buforowanie na potrzeby wycinania fragmetow danych
s_wynik0=size(wynik0);

%5a) Sprawdzenie warunku brzegowego z lewej strony (czy poczatek sygnalu s nie
%przypada w srodku spike'a)
if s_wynik0(2)>0
    for j=1:min(10,s_wynik0(2))            
        if wynik0(1,j)<left
            wynik=wynik0(:,j+1:s_wynik0(2));
        end
    end
end
wynik0=wynik; %buforowanie na potrzeby wycinania fragmetow danych
s_wynik0=size(wynik0);

%5b) Sprawdzenie warunku brzegowego z prawej strony (czy koniec sygnalu s nie
%przypada w srodku spike'a)
if s_wynik0(2)>0
    for j=s_wynik0(2):-1:max(s_wynik0(2)-10,1)            
        if wynik0(1,j)>dl-right
            wynik=wynik0(:,1:s_wynik0(2)-j-1);
        end
    end
end
    
clear wynik0;
s_wynik=size(wynik);
    
blad=zeros(3,s_wynik(2));
en=zeros(1,s_wynik(2));
ilosc=s_wynik(2);

%6. Wyznaczanie parametrow bledu dla kolejnych spike'ow
for j=1:s_wynik(2);
    spike=[wynik(1,j)-left:wynik(1,j)+right];
    blad(1,j)=min(roznica(spike));
    blad(2,j)=max(roznica(spike));
	ofset=mean(s(spike));
	ofset=0;
	blad(3,j)=energia(roznica(spike),ofset);
	en(1,j)=energia(s(spike),ofset);
end

%size(blad)
%figure(figura);
%plot(sqrt(en(1,:)),sqrt(blad(3,:)),'bd');
%axis([0 600 0 20]);
y=[wynik(1,:)' sqrt(en(1,:))' sqrt(blad(3,:))']';
%end