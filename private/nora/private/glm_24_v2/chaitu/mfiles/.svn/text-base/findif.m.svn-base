function fderr = findif(x, pars)
% call:  fderr = findif(x, pars)
% check function, gradient calculations provided by pars.fgname, pars.matfun
if size(x,2) ~= 1
   error('x must have only one column')
end
[f,g] = feval(pars.fgname, x, pars);
if f == inf
   error('findif: inital f is inf')
end
if size(g) ~= size(x)
   disp('size of gradient is wrong')
   keyboard
end
for k=1:16
   h = 10^(-k);
   d = rand(size(x));
   xpert = x + h*d;
   fpert = feval(pars.fgname, xpert, pars);
   fdif(k) = (fpert-f)/h;
   gd = g'*d;
   fderr(k) = abs(fdif(k) - gd);
end
figure(1)
semilogy(fderr,'*')
