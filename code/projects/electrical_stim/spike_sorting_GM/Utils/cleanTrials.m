function data=cleanTrials(data)
% cleanTrials erases the first trace of each condition j
% to avoid undesirable weird traces
% input:    data{j,e}
% output:   data{j,e} only taking the rows from 2 to I(j)
% Gonzalo Mena 06/2015

J = size(data,1);
E = size(data,2);
for j = 1:J
    for e=1:E
    dataAux{j,e} = data{j,e}(2:end,:);
    end
end
data=dataAux;
