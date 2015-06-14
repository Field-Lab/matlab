function dataVecJ=collapseDataE(data)
% collapseDataE concatenates the data, originally as data{j,i}
% as a J dimensional cell array dataVecJ{j}, the data is
% collapsed horizontally
%   input: data{j,e}
%   output: dataVecJ{j};
%   Gonzalo Mena 06/2015
J = size(data,1);
E = size(data,2);

for j = 1:J
    
    dataVecJ{j}=[];
    
    for e= 1:E
        dataVecJ{j}=[dataVecJ{j} data{j,e}];
    end
end