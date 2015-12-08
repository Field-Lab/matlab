function[rho theta] = normalize_maxspikes(rho, theta)

%Function that normalizes by the direction with the maximum spike rate
%Inputs are rho (the values) and theta(the angles)
%Function returns the normalized radius and angle matrix after removing the rows that divide
%by zeros and have Nans and Inf

%Sneha Ravi 
  %Last revision: 12-18-2012
  
for j = 1:length(rho)
    norm = [];
    norm = max(rho{j,1}')';
    norm = repmat(norm, 1, size(rho{j,1},2));
    rho{j,1} = rho{j,1}./norm;
%     [rho{j,1} theta{j,1}] = exciseRows(rho{j,1}, theta{j,1});
    rho{j,1} = exciseRows(rho{j,1});

end
end
