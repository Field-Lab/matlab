function wynik=detect2_fb(read_param,detect_param,block);
%wynik=detect2_fb(read_param,detect_param,block);
%Detekcja pikow dla duzych rekordow danych, z dodatkowym okreleniem
%wysokosci pikow. Dzialanie polega na 
%cyklicznym okreslaniu parametrow "read_param" i wywolywaniu
%funkcji '-> detect2_f' z tym parametrem. Argument "block" okresla
%maksymalny rozmiar bloku danych, na ktorych procedura "detect_f"
%bedzie dokonywala obliczen (zalecane takie okreslenie wartosci
%zmiennej "block", dla ktorej uniknie sie korzytania z pamieci
%wirtualnej).
%Wyniki zwracane przez funkcje moga w niewielkim stopniu zalezec
%od wartosci argumentu "block" ze wzgledu na ignorowanie przez
%funkcje 'detect2_f' pikow znajdujacych sie na granicach przedzialu.

samples=read_param.samples;
dlugosc=samples(2)-samples(1)+1;

ilosc=floor(dlugosc/block);
reszta=dlugosc-ilosc*block;

wynik=zeros(3,0);

for i=1:ilosc
    smps=[samples(1)+(i-1)*block samples(1)+i*block-1]
    read_param.samples=smps;
    w=detect2_f(read_param,detect_param);
    ws=size(w);
    if ws(2)>1
      w(1,:)=w(1,:)+(i-1)*block;
      w(2,:)=w(2,:)+(i-1)*block;
      wynik=[wynik w];
    end
end

if reszta>0
    smps=[samples(1)+ilosc*block samples(2)]   % ???
    read_param.samples=smps;
    w=detect2_f(read_param,detect_param);
    w(1,:)=w(1,:)+ilosc*block;
    w(2,:)=w(2,:)+ilosc*block;
    wynik=[wynik w];
end