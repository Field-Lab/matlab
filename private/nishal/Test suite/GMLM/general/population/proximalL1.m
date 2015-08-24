function xnew = proximalL1(x,lambda)
% Do soft thresholding

xnew = 0*x;
xnew(x<=-lambda) = x(x<= -lambda) + lambda;
xnew(x>=lambda) = x(x>=lambda) - lambda;

end