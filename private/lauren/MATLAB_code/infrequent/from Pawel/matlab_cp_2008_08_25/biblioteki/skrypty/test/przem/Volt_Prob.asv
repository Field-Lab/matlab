clear

cd H:/pliki/nauka/stymulacja/chip/testy/2006_03_24/

dane1=dlmread('awg',',');        %dane z oscyloskopu
dane2=importdata('awgp');          %parametry pracy - skala napiecia i czasu

volty=((4*dane2(1)/32000))*dane1;   %wyniki pomiarow w Voltach

maksimum=max(volty);                %maksymalna wartoœæ, zmienna pomocnicza

N=max(size(volty));                 %N - rozmiar tablicy danych volty
k=1                                 %k-zmienna pomocnicza w iteracjach

for i=1:N                           %-petla 
    if(volty(i)>(0.5*maksimum))
        k=k+1;
    end     
end
 
peak=zeros(1,k);
k=1;

for i=1:N
    if(volty(i)>(0.5*maksimum))
    peak(k)=volty(i);
    k=k+1;
    end
end


 srednia=mean(peak);
 plot(peak);
 k;
 srednia;
end