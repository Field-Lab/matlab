function data=cleanTrials(data)

J = size(data,1);
E = size(data,2);
for j = 1:J
    for e=1:E
    dataAux{j,e} = data{j,e}(2:end,:);
    end
end
data=dataAux;
