function s=readcnst(read_param);
%Funkcja czytajaca dane z pliku (JEDEN KANAL), z odjeciem trendu
%dlugoczasowego.
%Parametr wejsciowy "read_param" jest struktura o nastepujacych 
%polach:
%1. name - nazwa pliku. W katalogu biezacym musi sie znajdowac plik
%z danymi oraz plik z trendem dlugoczasowym. Pierwsze 7 znakow 
%w nazwach obu plikow musi byc identyczne i zgodne z argumentem 
%"name". Plik z danymi zawiera ponadato przyrostek "conv"
%(nie rozszerzenie - tzn. bez kropki), a plik z trendem - przyrostek
%"cnst";
%2. header - dlugosc naglowka w pliku "*conv" (w bajtach);
%3. nrchns - ilosc kanalow (w plikach "*conv" i "*cnst" identyczna);
%4. channels - numer kanalu, dla ktorego dokonuje sie odczytu;
%5. samples - wektor dwuelemetowy (jednowierszowy), okreslajacy
%pierwsza i ostatnia probke do zaladowania.

name=read_param.name;
header=read_param.header;
nrchns=read_param.nrchns;
channels=read_param.channels;
samples=read_param.samples;

lchannels=length(channels);

name0(1,1:7)=name(1,1:7);
name0(1,8:11)='conv';

name1(1,1:7)=name(1,1:7);
name1(1,8:11)='cnst';

fid0=fopen(name1,'r');
window=fread(fid0,1,'int32');
marg=fread(fid0,1,'int32');

if window<1
    error('Dlugosc okna w pliku .cnst musi byc dodatnia');
end

if marg<0
    error('Dlugosc marginesu w pliku .cnst nie moze byc ujemna');
end

fclose(fid0);

l_dane=raportconv(name0,header,nrchns); %dlugosc kanalu w *conv

if samples(1)<1
    error('Zbyt male wartosci w "samples"');
end

if samples(2)>l_dane
    error('Zbyt duze wartosci w "samples"');
end

ilosc=l_dane/window;

if round(ilosc)~=ilosc
    warning('Dlugosc bloku danych nie jest calkowita wielokrotnoscia dlugosci okna');
end

ilosc=floor(ilosc);

l_cnst=raportconv(name1,12,nrchns); %ilosc wartosci trendu w plik .cnst;

offset1=floor((ilosc-l_cnst)/2);  %ile wartosci trzeba dopisac na poczatku wektora
                                  %wartosci trendu;
offset2=ilosc-offset1;            %j.w. ale na koncu wektora;

y0=readconv(name1,12,nrchns,channels,[1 l_cnst]);

trend(1,(offset1+1):(offset1+l_cnst))=y0;
trend(1,1:offset1)=y0(1,1);
trend(1,(offset1+l_cnst+1):ilosc)=y0(1,l_cnst);  %uzupelnienie tablicy trendu na poczatku i koncu;

s=readconv(name0,header,nrchns,channels,samples);

start=floor((samples(1)-1)/window);  %ktore wartosci z pliku .cnst beda uzywane;
stop=floor((samples(2)-1)/window);

for i=start:stop
    start1=max(i*window+1,samples(1))-samples(1)+1;
    stop1=min((i+1)*window,samples(2))-samples(1)+1;
    s(1,start1:stop1)=s(1,start1:stop1)-trend(i+1);
end
    

















