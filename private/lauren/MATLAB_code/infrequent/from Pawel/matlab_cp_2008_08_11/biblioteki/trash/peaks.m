function w=peaks(s,offset,prog,histereza,opoznienie);

s=s-mean(s);
%s=delconst(s,1000,4000);
y=detekcja2(s,0,prog,histereza);

l=length(y)
w0=zeros(l,3);
w0(:,1)=y';

for i=1:l
    s0=s(1,y(i):(y(i)+opoznienie))-(prog-histereza);    
    %size(s1);
    s1=diff(sign(s0));
    n=min(min(find([s1 3])));
    if n==opoznienie+1
        i
        error('Zbyt mala wartosc opoznienia');
    end
    
    w0(i,2)=n;
    w0(i,3)=max(s0(1,1:n));
end

w=w0;