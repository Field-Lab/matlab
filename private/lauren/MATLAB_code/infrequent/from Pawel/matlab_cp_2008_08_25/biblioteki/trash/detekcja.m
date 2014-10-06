function p=detekcja(x,offset,prog,histereza);
%Dokonuje detekcji sygnalow w rekordzie 'x',
%Tablica wynikowa zawiera trzy wiersze - 
%kolejno pocztek piku (nr probki, dla ktorej wykryto 
%przekroczenie progu), jego szerokosc (az do zejscia ponizej progu)
%i wysokosc wzgledem wartosci offset.

a=size(x);
numer=0;
flaga=0;

x=x-offset;

piki=zeros(3,a(2));

for j=1:a(2)
    if flaga==0
        if abs(x(1,j))>(prog+histereza(2))
            flaga=1;
            numer=numer+1;
            piki(1,numer)=j;
        end
    end

    if (flaga==1)
        if abs(x(1,j))<(prog+histereza(1))
            flaga=0;
            piki(2,numer)=j;
            if(sign(x(1,piki(1,numer)))==1)
                piki(3,numer)=max(x(1,piki(1,numer):piki(2,numer)));
            else
                piki(3,numer)=min(x(1,piki(1,numer):piki(2,numer)));
            end
        end
    end
    
end
    
ilosc=numer;
p=piki(:,1:numer);