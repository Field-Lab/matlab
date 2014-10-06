function y=phist(s,n);

a=rand(1,length(s))-0.5;
if n~=0 
    hist(s+a,n);
    y='ok';
else
    hist(s+a);
    y='ok';
end

