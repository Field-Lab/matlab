function r=inveyefilter(n,f3,f1,f2,fp,f0);
%wyznacza odpowiedz impulsowa filtru odwrotnego do filtru
%pasmowo-przepustowegotypu rc2cr2. 
%Dane wejsciowe: N - rzad filtru, f3 - czest. gran. odsprzezenia 
%stalopradowego, f1 - czest.
%graniczna sekcji gornoprzep., f2 - czest. graniczna sekcji
%dolnoprzepustowej (f1<f2), fp - czestotliwosc probkowania,
%f0 - czest., ponizej ktorej charakterystyka jest nieistotna
%(dla unikniecia osobliwosci przy f=0). 
%Wyjscie: odpowiedz impulsowa zawierajaca 2N+1 probek
%(a wiec dlugosc zawsze nieparzysta).

N=2*n+1;
df=fp/N;
k=ceil(f0/df)     %ponizej f0 charakterystyka jest przyblizana
f=[0:(N-1)]/N*fp;
ch_cz=zeros(1,N);
%df
omega=2*pi*f;
t1=1/(2*pi*f1);
t2=1/(2*pi*f2);

j=sqrt(-1); %dla formalnosci
%j=1;
ch_cz=(j*omega*t1)./((1+j*omega*t1).*(1+j*omega*t2));
ch_cz=ch_cz.*ch_cz; %bo drugiego rzedu;

%Mamy filtr, teraz stala na wejsciu:
t3=1/(2*pi*f3);
ch_cz3=(j*omega*t3)./(1+j*omega*t3);
ch_cz=ch_cz.*ch_cz3;

for i=1:k 
   ch_cz(1,i)=ch_cz(k+1); %przyblizanie dla f<f0
end
ch_cz=ch_cz.^(-1); %bo odwrotny
ch_cz=ch_cz.*exp(j*omega*n/fp); %korekta fazy
for i=(n+2):N
   ch_cz(1,i)=conj(ch_cz(1,N+2-i));
end
%r=real(ifft(ch_cz));
r=real(ifft(ch_cz));
okno=blackman(2*n+1)';
r=r.*okno;

