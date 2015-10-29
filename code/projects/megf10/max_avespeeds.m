function[magave] = max_avespeeds(mag)

%Function returns average of magnitude vectors calculated across all temporal periods

% Sneha Ravi 
% Last revision: 12-18-2012

for i = 1:length(mag)
    magave(i,1:length(mag{i,1})) = mag{i,1};
end
magave = mean(magave);
end