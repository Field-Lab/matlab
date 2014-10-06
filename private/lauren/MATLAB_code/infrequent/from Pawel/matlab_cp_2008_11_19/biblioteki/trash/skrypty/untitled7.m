a=rand(1,1000000);
b=rand(1,1000000);
przedzial=[50000 50000];

wykl=[12:1:20];

for i=1:length(wykl)
    if (2^wykl(i)>(przedzial(1)+przedzial(2)))
        clear c;
        wykladnik=wykl(i)
        s=clock;
        s(6);
        d=fastcorr2(a,b,przedzial,wykl(i));
        s1=clock;
        s1(6)-s(6)
    end
end

'koniec'