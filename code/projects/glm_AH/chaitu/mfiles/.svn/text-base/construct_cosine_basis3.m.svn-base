% functional form of raised cosine (with logarithmic time scale) is
% g_j(t) = (1/2)cos(a*log[b*t + c] - phi_j) + (1/2)
% for t such that (phi_j - pi) <= a*log(b*t + c) <= (phi_j + pi)

% alpha - starting point (can be negative)
% beta - finishing point
% b - time scaling INSIDE log
% nfilters - number of filters to construct
% dt - timestep (frame or spike -- note this parameter is working weird..)

% based on construct_cosine_basis2.m, edoi, 2011-01-07.
% add: in the first basis, basis element is set 1 until the peak (as in Pillow et al, 2008, suppl.)

function [basis, phi,tvec] = construct_cosine_basis3(alpha,beta,b,nfilters,dt,spacing,fl_out)
if ~exist('fl_out','var')
   fl_out = 0;
end

tvec = alpha:dt:beta;
% define c s/t [b*t + c] = dt when t = alpha; this is necessary for
% the logarithmic time scaling (where t need to be > 0)
c = -b*alpha + dt;
if (~exist('spacing','var'))
   spacing = pi/2; % default spacing (square sum is constant)
end
a = (spacing*(nfilters-1) + 2*pi) * 1/log((b*beta + c)/(b*alpha + c));
minphi = a*log(b.*alpha + c) + pi;
maxphi = a*log(b.*beta  + c) - pi;

if (maxphi < minphi)
   fprintf('ERROR: alpha beta are unsuitable for cosine basis construction, minphi=%f, maxphi=%f !\n',...
      minphi,maxphi);
   return;
end

T = length(tvec);

phi = linspace(minphi,maxphi,nfilters);
if (length(phi) > 1)
   A = norm(diff(phi) - repmat(spacing,1,nfilters-1));
   if (A > 10^(-6))
      fprintf('ERROR: spacing is incorrect to precision %f!n',A);
   end
end

basis = zeros(T,nfilters);
% construct the postspike basis functions
arg = a.*log(b.*tvec + c);
for k=1:nfilters
   idx_k =  find((arg > phi(k)-pi) & (arg < phi(k)+pi));
   if isempty(idx_k)
      fprintf('timespan of filter %d is empty!\n',k);
   else
      if fl_out == 1
         fprintf('timespan of filter %d is (%f,%f)\n',k,tvec(idx_k(1)),tvec(idx_k(end)));
      end
      basis(idx_k,k) = 0.5.*(cos(arg(idx_k) - phi(k))+1);
      if k == 1 % add
         idx_k0 =  (arg < phi(1));
         basis(idx_k0,1) = 1;
      end
   end
end

% correction of the initial sharp raise: edoi, 2012-01-20
nbin_ar = find(tvec<3/4*10^-3,1,'last'); % for modeling the absolute refractory period
if sum(idx_k0) > nbin_ar
   nbin_shift = sum(idx_k0) - nbin_ar;
   basis = [basis(nbin_shift+1:end,:);zeros(nbin_shift,nfilters)];
end



