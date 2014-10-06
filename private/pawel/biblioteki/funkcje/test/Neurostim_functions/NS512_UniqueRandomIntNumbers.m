function Numbers=NS512_UniqueRandomIntNumbers(Range,ArrayLength);

Numbers=zeros(1,ArrayLength);
Values0=[1:Range];
Values=Values0;
for i=1:ArrayLength
    %Values=[1:Range]
    ValueID=ceil(rand(1,1)*Range);
    Numbers(i)=Values(ValueID);
    
    Values(ValueID)=0;
    Values=nonzeros(Values);
    Range=Range-1;
end