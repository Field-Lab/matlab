function [r,m,s] = average_radius(datarun, cell_specification)
%
% greschner


% get cell numbers
index=get_cell_indices(datarun,cell_specification);


if length(index)<10
    warning(sprintf('average_radius: mean sd computed on small numbers: %d', length(index)))
end
    
temp=zeros(size(index));
for i=1:length(index)
    temp(i)=sqrt(prod(datarun.(datarun.default_sta_fits).sta_fits{index(i)}.sd));
end

r=robust_mean(temp);
m=mean(temp);
s=std(temp);

















