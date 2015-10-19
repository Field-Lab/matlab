function[dsindex] = DS_index_one(mag, null)

% Function calculates and returns direction selective index, dsi = Vector-Null/Vector+Null

% Sneha Ravi 
% Last revision: 12-18-2012

dsindex = cell(length(mag),1);
for i = 1:length(mag)
    dsindex{i,1} = (mag{i,1}-null{i,1})./(mag{i,1}+null{i,1});
%     [dsindex{i,1} ntimp] = exciseRows(dsindex{i,1}', zeros(size(dsindex{i,1},2),1));
    dsindex{i,1} = exciseRows(dsindex{i,1}');
end
end