function [slope,lin_err,dane]=dac_lin3(dac,dane0);
dane=[dane0(1,1:127) (dane0(1,128)+dane0(1,129))/2 dane0(1,130:256)];

[p,s] = POLYFIT(dac,dane,1);
y = POLYVAL(p,dac);

aa=[dane./dac];                 %aa - tablica wspolczynnikow nachylenia dla wszystkich wartosci
aa1=[aa(1,1:127) aa(1,129:255)];
a1=min(abs(aa1))               %a1 - najmniejszy wspolczynnik nachylenia
a2=max(abs(aa1))               %a2 - najwiekszy wspolczynnik nachylenia

N=500;
step=(a2-a1)/N;

glob_emax=10;                   %glob_emax - najwieksza wartosc bledu maksymelnego dla najlepiej dopasowanego wspolczynnika nachylenia

for i=1:N
    c=a1+i*step; 
    y=c*dac;        
    err=abs(y-dane);    
    emax=max(err);              %emax - najwiekszy blad dla danego wspolczynnika nachylenia
    if (emax < glob_emax)
        glob_emax=emax;
        a_bestfit=i;            %a_bestfit - numer najlepiej dopasowanego wsp. nachylenia
    end
    
 end 
 
c=a1+a_bestfit*step;

slope=c;
y=slope*dac;
lin_err=(y-dane)./slope;