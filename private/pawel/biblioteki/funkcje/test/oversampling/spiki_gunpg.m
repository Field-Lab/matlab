function [p,filenames,spikes,detect_param]=spiki_malpa(type,ktore,system);
%Funkcja zwraca namiary na dobre spiki dla guinepig. Podaje sie typ
%spikow (1-cell body, inne - axonal), funkcja zwraca: sciezke do odpowiednich plikoe (wszystkie pliki z
%maja 2004), nazwy plikow, strukture zawierajaca - dla kazdego spika -
%nazwe pliku i wspolrzedne spika, oraz parametry detekcji uzywane przy
%wyszukiwaniu spikow.
%system - linux=1, windows=2 (inna sceizka dostepu!)

%nie zmieniac parametrow detekcji!!!! numery spikow sa prawidlowe TYLKO dla
%prog=60, histereza=30, znak=-1.
'fge'

prog=60; 
histereza=30;
znak=-1;
detect_param=struct('prog',prog,'histereza',histereza,'znak',znak);
 
if system==1
	p='/mnt/data4/dane/2004/26maja/2003-12-08-0-001'
else
	p='K:\dane\2004\26maja\2003-12-08-0-001'
end

cd(p);
filenames={'electrode45.txt' 'electrode137.txt' 'electrode148.txt' 'electrode155.txt' 'electrode169.txt' 'electrode221.txt' 'electrode298.txt' 'electrode362.txt' 'electrode416.txt' 'electrode499.txt'};
nr_of_files=length(filenames);


%1 spiki typu "cell body"
cb_spikes=zeros(19,2);
cb_spikes(:,1)=[1 2 2 3 4 4 5 5 6 6 7 7 8 8 9 9 9 10 10]';
%cb_spikes(:,2)=[85 30 240 227 2 63 158 181 10 147 1 20 4 5 61 68 167 14 17]';
cb_spikes(:,2)=[1132 30 240 227 2 63 158 181 10 147 1 20 4 5 61 68 167 14 17]';

%2. Spiki typu "axonal"
ax_spikes(:,1)=[1 1 2 2 3 3 4 5 5 7 7 8 8 8 8 9 9 9]';
ax_spikes(:,2)=[185 300 231 234 106 254 285 5 285 227 228 3 225 227 281 80 161 180]';


if type==1
    spikes=cb_spikes;
else
    spikes=ax_spikes;
end

spikes=spikes(ktore,:);

for j=1:nr_of_files %dla kazdego pliku:
    if find(spikes(:,1)==j)'	    
        filenames{j}
        s=importdata(filenames{j})';
        s=s-mean(s);
    
        y=detect2(s,detect_param);    
        positions=find(spikes(:,1)==j)' %numery dobrych spikow w tym pliku    
        coordinates(positions)=y(1,spikes(positions,2)); %wpisanie wspolrzednych
    end	
end

spikes(:,2)=coordinates';

%el_10: cb - 35,66 axonal - 64,197,207
%el167: cb - 1,25,141,281,242; axonal - 222,229,9,10,11,12,13,264
%el168: cb - 113,138; axonal - 1,2,285
%el190: cb - 244,290,295,20,80,162 (duuzo, b.wysokie); axonal - 126,129; dziwne -191,192,287,288
%el200: cb - 202;axonal - 8,26,89,92,224,244; dziwne: 54
%el353: cb - 196,197; axonal - 17,  dziwne:201,63,247,55
%el399: cb - 192,252,292 axonal - 105,186; dziwne - 26,68,151,248,51;
