function w=shwpks(s,peaks,samples,numbers,fig);
%w=shwpks(s,peaks,samples,numbers,fig);
%Funkcja rysujaca okreslone piki. 
%s - sygnal (wektor);
%peaks - wartosci poczatkowe pikow (przekroczenie progu);
%samples - dwa elementy, okreslajace ile probek przed oraz po
%          poczatku piku ma byc wyswietlane (np. [20 50]) - 20
%          na lewo, 50 na prawo - w sumie 71;
%numbers - numery pikow (wg kolejnosci w wektorze 'peaks';
%fig - numer 'figury', na ktorej funkcja bedzie rysowac;

figure(fig);
clf(fig);
l=length(numbers);
p=sqrt(l);

a2=floor(p);
%a2=20;

if a2==0
   a2=1;
end

a1=ceil(l/a2);

figure(fig);

for i=1:l
    subplot(a1,a2,i);
    t=[(peaks(1,i)-samples(1)):(peaks(1,i)+samples(2))];
    s0=s(1,t(1,1):t(1,length(t)));
    %plot(t,s0);
    plot(t,s0);
    grid on;
end

w=l;