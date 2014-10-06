function [p,filenames,spikes,detect_param]=spiki_gnpg(type);
%Funkcja zwraca namiary na dobre spiki dla guinepig. Podaje sie typ
%spikow (1-cell body, inne - axonal), funkcja zwraca: sciezke do odpowiednich plikoe (wszystkie pliki z
%maja 2004), nazwy plikow, strukture zawierajaca - dla kazdego spika -
%nazwe pliku i wspolrzedne spika, oraz parametry detekcji uzywane przy
%wyszukiwaniu spikow.

%nie zmieniac parametrow detekcji!!!! numery spikow sa prawidlowe TYLKO dla
%prog=60, histereza=30, znak=-1.
prog=60; 
histereza=30;
znak=-1;
detect_param=struct('prog',prog,'histereza',histereza,'znak',znak);
 
p='I:\dane\2004\26maja\2003-12-08-0-001';
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

for j=1:nr_of_files %dla kazdego pliku:

    s=importdata(filenames{j})';
    s=s-mean(s);
    
    y=detect2(s,detect_param);    
    positions=find(spikes(:,1)==j)' %numery dobrych spikow w tym pliku    
    coordinates(positions)=y(1,spikes(positions,2)); %wpisanie wspolrzednych
end

spikes(:,2)=coordinates';

%1. Spiki typu "cell body":
%I:\dane\2004\26maja\2003-12-08-0-001 - kolejne pliki:
%electrode45: cell body: 85, axonal: 185,300
%electrode137: cell body: 30, 240 axonal:231, 234
%electrode148: cell body 227, axonal: 106, 254
%electrode155: 2,63 axonal: 285
%electrode169: cb: 158, 181 axonal: 5,285
%el221: cb - 10,147 axonal:none
%el298: cb - 1,20; axonal:227, 228
%el362: cb - 4,5; axonal: 3, 225, 227, 281 (mnooostwo!)
%el416: cb - 61,68,167; axonal: 80, 161, 180 (duzo)
%el499: cb - 14, 17; axonal:none
