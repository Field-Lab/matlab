function y=delconst(x,window,marg);

l=length(x);
n=floor(l-2*marg)/window
y=zeros(1,l-2*marg);

for i=1:n
    i
    start=(i-1)*window+marg+1;
    stop=i*window+marg;
    a=mean(x(1,(start-marg):(start+marg)));
    size(a);
    y(1,((i-1)*window+1):i*window)=x(1,start:stop)-a;
end
