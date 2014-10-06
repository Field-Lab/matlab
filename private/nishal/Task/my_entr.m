function y=my_entr(x)
y=[];
for i =1:size(x,1)
    
    if(x(i)==0 || x(i)==1 )
        y(i)=0;
    else
        y(i)= -x(i).*log2(x(i));
    end
end

y=y';
end