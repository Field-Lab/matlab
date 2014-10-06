N=100;
a=zeros(1,N);
s=0;

for i=1:N
    a(1,i)=N/(N-i+1);
    s=s+a(1,i);
end