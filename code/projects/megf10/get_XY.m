function[X Y] = get_XY(rho, theta)

%Function converts magnitude of spike numbers in each direction to X and Y polar coordinates


%Sneha Ravi 
  %Last revision: 12-18-2012

X = cell(length(rho), 1);
Y = cell(length(rho), 1);

for i = 1:length(rho)
    [X{i,1} Y{i,1}] = pol2cart(theta{i,1}, rho{i,1}); 
end
end