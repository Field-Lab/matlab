function [x, it, err] = newton(f, df, x0, tol, nmax)
  % compute a root of f using Newton's method, starting at x0

  % default values for some arguments
  if nargin < 4
    tol = 1e-5;
  end
  if nargin < 5
    nmax = 40;
  end

  x = x0;
  for it = 1:nmax
    s = df(x) \ f(x);
    x = x - s;

    % stop if we are below the tolerance
    if (norm(s) <= tol*norm(x))
      return;
    end
  end
end