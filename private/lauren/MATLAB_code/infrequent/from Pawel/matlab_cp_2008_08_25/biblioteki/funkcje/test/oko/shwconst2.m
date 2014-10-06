function y=shwconst2(x,window,marg,prog);

l=length(x);
n=floor(l-2*marg)/window;
y=zeros(1,n);

for i=1:n
    i;
    start=(i-1)*window+marg+1;
    stop=i*window+marg;
    s=x(1,(start-marg):(start+marg));
    a=mean(s);
    f=find(abs(s-a)>prog);
    
    s(1,[f])=a;
    a=mean(s);
    
    %f length(f)==0
        %r=0;
        %else
        %r=sum(s([f])-a)/length(f);
        %end
    %y(1,i)=a-sum(s([f])-a)/length(f);
    y(1,i)=a;
    %y(1,((i-1)*window+1):i*window)=s(1,(marg+1):(marg+window))-a;
end
