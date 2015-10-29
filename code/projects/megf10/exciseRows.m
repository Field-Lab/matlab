
function [rho] = exciseRows(rho)
%Function removes any rows that contain Inf or Nans from rho matrix
%Input: Rho Matrix that contains Inf and Nans
%Output: Rho Matrix without Inf or Nans as well as adjusted Theta Matrix

%Sneha Ravi 
  %Last revision: 12-18-2012


%l1 = size(rho,1);
rho(any(isnan(rho),2),:) = 0;
rho(any(isinf(rho),2),:) = 0;
%l2 = size(rho,1);
%theta(l2+1:l1,:) = [];
end