function x = generate_ica_dist(type, points)
% GENERATE_ICA_DIST    Generates distribution of data for ICADEMO
%

switch type
 case 1
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % exponential distribution
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % uniform distribution
  x = rand(1,points);

  % binary distribution
  bin = sign(rand(1,points) - 0.5);

  % inverted exponential function
  % Note: zero almost never happens
  alpha = 8;
  x = bin .* 1/alpha .* log(1-x);

 case 2
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % gaussian distribution
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % standard deviation of 0.2
  x = randn(1,points) * 0.2;

 case 3
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % uniform distribution
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  x = rand(1,points) - 0.5;

  % squish it to screen size
  x = x * 0.8;

end


