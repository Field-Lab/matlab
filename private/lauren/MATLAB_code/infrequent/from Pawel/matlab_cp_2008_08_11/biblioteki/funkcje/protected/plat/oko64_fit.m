function parametry=oko64_fit(dane,funkcja,start,krok,prog);

d_size=size(dane);
parametry=zeros(d_size(1),length(start));

t=[1:d_size(2)];
for i=1:d_size(1)
%for i=1:10
  s=dane(i,:);
  i;
  max(s);
  min(s);
  prog;
  if (max(s)-min(s))>prog
    i;
    
    %
    %f=fft(s-mean(s));
    %start(1,3)=angle(f(2));
    %start(1,1)=max(max(abs(f)))/length(s)*2;
    %start(1,2)=
    %
    start=sinus_fit_pocz(s);
 
    [b,p]=nielfit2(t,s,'sinus_fit',start,krok);
    
    bsize=size(b);
    b;
    psize=size(p);
    i;
    parametry(i,:)=b;
  end
end
