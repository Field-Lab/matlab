function wynik=raport(nazwa);
%clear;

%nazwa='probka3_filt.bin';
fid1 = fopen(nazwa,'r');

fp=fread(fid1,1,'int')
ilosc=fread(fid1,1,'int')
if (ilosc>64)
   error('Ilosc kanalow wieksza od 64');
end
kanaly=fread(fid1,ilosc,'int')
dlugosc=fread(fid1,1,'int')
offset=fread(fid1,32,'int')

fclose(fid1);

wynik='ok.';