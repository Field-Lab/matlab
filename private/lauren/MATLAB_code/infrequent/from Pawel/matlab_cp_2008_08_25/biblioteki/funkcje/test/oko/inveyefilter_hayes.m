function [filter,inv_filter]=inveyefilter_hayes(n,f0,f1,f2,fp);
%wyznacza odpowiedz impulsowa filtru odwrotnego do filtru
%pasmowo-przepustowegotypu rc2cr2. 
%Dane wejsciowe: N - rzad filtru, f0 - czest. gran. odsprzezenia 
%stalopradowego, f1 - czest.
%graniczna sekcji gornoprzep., f2 - czest. graniczna sekcji
%dolnoprzepustowej (f1<f2).
%Wyjscie: odpowiedz impulsowa zawierajaca 2N+1 probek
%(a wiec dlugosc zawsze nieparzysta).

N=2*n+1;
df=fp/N;

f=[0:(N-1)]/N*fp;
ch_cz=zeros(1,N);
%df
omega=2*pi*f;
%stale czasowe:
t0=1/(2*pi*f0);
t1=1/(2*pi*f1);
t2=1/(2*pi*f2);

j=sqrt(-1); %dla formalnosci
%j=1;
ch_cz=(j*omega*t1)./((1+j*omega*t1).*(1+j*omega*t2));
ch_cz=ch_cz.*ch_cz; %bo drugiego rzedu;

%Mamy filtr, teraz stala na wejsciu:
ch_cz0=(j*omega*t0)./(1+j*omega*t0);
ch_cz=ch_cz.*ch_cz0;

ch_cz=ch_cz.*exp(j*omega*n/fp); %korekta fazy
for i=(n+2):N
   ch_cz(1,i)=conj(ch_cz(1,N+2-i));
end

%filter=ifft(ch_cz);
filter=real(ifft(ch_cz));
[inv_filter,err]=spike(filter,N-1,N);
inv_filter=inv_filter';
%r=r.*okno;
