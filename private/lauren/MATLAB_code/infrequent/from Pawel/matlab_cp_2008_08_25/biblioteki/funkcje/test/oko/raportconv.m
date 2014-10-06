function dlugosc=raportconv(name,header,nrchns);

fid0=fopen(name,'r');

a=fseek(fid0,0,1); %skok na koniec pliku
p=ftell(fid0);

n=floor((p-header)/(2*nrchns));

if p>(n*2*nrchns+header)
    warning('Jezeli parametry wejsciowe funkcji sa prawidlowe, kanaly nie sa jednakowej dlugosci');
end

dlugosc=n;