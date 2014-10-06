function y=shwconst(name,channels,window,marg,len,header,nrchns);
%Funkcja wyznaczjaca przebieg dryfu w poszczegolnych 
%kanalach. 
%name - nazwa pliku;
%channels - numery kanalow;
%window - dlugosc okna czasowego (ilosc probek), od ktorego
%         odejmuje sie srednia w kazdej iteracji;
%marg - ilosc probek na lewo i prawo (jedna liczba) od okna,
%       ktore uwzglednia sie przy liczeniu sredniej; laczna 
%       ilosc probek aynosi window+2*marg;
%numbers - numery (poczatkowy i koncowy) okien czasowych. 
%          Nalezy pamietac, ze funkcja z definicji odrzuca 
%          na poczatku i koncu przebiegu ilosc probek 
%          rowna wartosci 'marg'.

l=length(channels);
ilosc=floor(len-2*marg)/window;
%n=floor(l-2*marg)/window
start=(numbers(1)-1)*window+1;
stop=numbers(2)*window;
y=zeros(l,

for i=1:l
    r=readconv(name,header,nrchns,channels(i),[1 len]);
    
    
    

for i=1:n
    i;
    start=(i-1)*window+marg+1;
    stop=i*window+marg;
    a=mean(x(1,(start-marg):(start+marg)));
    y(1,((i-1)*window+1):i*window)=x(1,start:stop)-a;
end
