function dataVecJ=collapseDataE(data)

J = size(data,1);
E = size(data,2);

for j = 1:J
    
    dataVecJ{j}=[];
    
    for e= 1:E
        dataVecJ{j}=[dataVecJ{j} data{j,e}];
    end
end