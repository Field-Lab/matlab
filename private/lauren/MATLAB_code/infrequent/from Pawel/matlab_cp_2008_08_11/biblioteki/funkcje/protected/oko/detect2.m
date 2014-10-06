function wynik=detect2_f(s0,detect_param);
%wynik=detect2_f(read_param,detect_param);
%
%Funkcja pokrewna do funkcji '-> detect_f', w stosunku do ktorej
%dodatkowo dokonywane jest okreslenie wysokosi kazdego piku. Funkcja
%zwraca tablice o rozmiarach 3xN, a odpowiednie elementy trzeciego
%wiersza okreslaja amplitude danego piku.
%Zobacz takze: detect_f, detect2_fb.

prog=detect_param.prog;
histereza=detect_param.histereza;
znak=detect_param.znak;

if (histereza>=prog)
    error('Histereza musi byc mniejsza od progu');
end

%s0=readcnst(read_param);

switch znak
case -1
    s=-s0;
case 0
    s=abs(s0);
otherwise
    s=s0;
end

dw=diff(s);
f1=find(s>(prog));
f1=[0 f1];

df1=diff(f1)-1;
f11=find(df1)+1;
w1=f1([f11]); %indeksy probek, dla ktorych nastepuje
              %przekroczenie wyzszego progu;
if length(w1)<1
    wynik=0;
    result='No spikes found'
    return;
end

f1=find(s>(prog-histereza));
f1=[f1 length(s)+1];

df1=diff(f1)-1;  
f11=find(df1);    %wskazania do okreslonych pikow
w2=f1([f11])+1;

w=sort([w1 w2]);
ds=dw([w]-1); %pochodne funkcji w punktach granicznych
              % (nachylenie zboczy przecinajacych progi)

sgn0=min(find(ds>0));
sgn1=max(find(ds<0));

w=w(1,sgn0:sgn1);
ds=ds(1,sgn0:sgn1);
sgn=sign([0 ds]);
dsgn=diff(sgn);
wynik(1,:)=w(find(dsgn>0));
wynik(2,:)=w(find(dsgn<0));
sw=size(wynik);

for i=1:sw(2)
    wynik(3,i)=max(s(1,wynik(1,i):wynik(2,i)))*sign(s0(wynik(1,i)));
end

ilosc=length(wynik);
