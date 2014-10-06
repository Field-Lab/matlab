function wynik=detekcja5(s,prog,histereza);
%Detekcja pikow dla sygnalu BEZ SKLADOWEJ STALEJ.
%s - signal;
%prog - prog dyskryminacji;
%histereza - o ile nizszy ma byc prog okreslajacy koniec piku.
%
%uwaga! Moze byc problem z ostatnim pikiem w sygnale.

s=abs(s);
dw=diff(s);
f1=find(s>(prog));
f1=[0 f1];

df1=diff(f1)-1;
%df1=[1 df1]
f11=find(df1)+1;
w1=f1([f11]); %indeksy probek, dla ktorych nastepuje
              %przekroczenie wyzszego progu;
if length(w1)<1
    wynik=0;
    return;
    %error('Nie wykryto pikow');
end

f1=find(s>(prog-histereza));
f1=[f1 length(s)+1];

df1=diff(f1)-1;  
%df1=[1 df1]
f11=find(df1);    %wskazania do okreslonych pikow
w2=f1([f11])+1;

w=sort([w1 w2]);
ds=dw([w]-1); %pochodne funkcji w punktach granicznych
              % (nachylenie zboczy przecinajacych progi)
%sgn=sign(ds);
sgn0=min(find(ds>0));
sgn1=max(find(ds<0));

w=w(1,sgn0:sgn1);
ds=ds(1,sgn0:sgn1);
sgn=sign([0 ds]);
dsgn=diff(sgn);
wynik(1,:)=w(find(dsgn>0));
wynik(2,:)=w(find(dsgn<0));

ilosc=length(wynik);
