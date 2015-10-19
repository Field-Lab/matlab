function[magmax] = max_allspeeds(mag)

%Function returns maximum of magnitude vectors over all temporal periods

% Sneha Ravi 
% Last revision: 12-18-2012

for i = 1:length(mag)
    magmax(i,1:length(mag{i,1})) = mag{i,1};
end
magmax = max(magmax);
end
