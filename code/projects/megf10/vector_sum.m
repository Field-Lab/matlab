function[U V angle mag] = vector_sum(X, Y)

%Function that sums spike numbers in all directions     

%Input: X and Y magnitudes of all cells for all temporal periods in all directions

%Output: Cartesian (U,V) and Polar (angle, mag) forms of calculated vector sum

% Sneha Ravi 
% Last revision: 12-18-2012


U = cell(length(X), 1);
V = cell(length(X), 1);
angle = cell(length(X), 1); 
mag = cell(length(X), 1); 

for i = 1:length(X)
    U{i,1} = sum(X{i,1}');
    V{i,1} = sum(Y{i,1}');
    [angle{i,1},mag{i,1}] = cart2pol(U{i,1}, V{i,1});  
end
end