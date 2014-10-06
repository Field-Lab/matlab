function Indexes=NS_IndexesForPairStimulation(number);

for i=1:7
    Indexes(i,i)=1;
    Indexes(i+7,i)=-1;
end

for i=2:7
    Indexes(14+(i-2)*8+1,i)=1;
    Indexes(14+(i-2)*8+1,1)=1;
    
    Indexes(14+(i-2)*8+2,i)=1;
    Indexes(14+(i-2)*8+2,1)=-1;
    
    Indexes(14+(i-2)*8+3,i)=-1;
    Indexes(14+(i-2)*8+3,1)=1;
    
    Indexes(14+(i-2)*8+4,i)=-1;
    Indexes(14+(i-2)*8+4,1)=-1;
    
    if i<7
        shift=1;
    else
        shift=-5;
    end
    
    Indexes(14+(i-2)*8+5,i)=1;
    Indexes(14+(i-2)*8+5,i+shift)=1;
    
    Indexes(14+(i-2)*8+6,i)=1;
    Indexes(14+(i-2)*8+6,i+shift)=-1;
    
    Indexes(14+(i-2)*8+7,i)=-1;
    Indexes(14+(i-2)*8+7,i+shift)=1;
    
    Indexes(14+(i-2)*8+8,i)=-1;
    Indexes(14+(i-2)*8+8,i+shift)=-1;
end