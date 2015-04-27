% Function that does preconditioned nonlinear CG. Polack-Ribiere and newton-raphson method are used
% for the linesearch

% Arguments:

% x0 - initial search point
% function handle that computes the function AND gradient at a point
% pars - parameter structure with the following fields:
%        - max_iter - max # of CG iterations
%        - tol - CG tolerance
%        - step - step size for the linesearch
%        - max_ls_iter - max # of linesearch iterations
%        - ls_tol - tolerance for linesearch
% basepars,stimpars,trainpars are the usual parameters for training the GLM
%
% precon_func - optional function that takes a vector and returns the
%               inverse of the preconditioner time that vector.

% Returns the optimal point xstar

function xstar = pcnlcg(x0,func,pars,basepars,stimpars,trainpars,precon_func)


if (~exist('precon_func','var'))
    precon_func = @(x) x;
end

% Initialize variables
i = 0;
k = 0;
x = x0;

% Initial function evaluation
[f g cifs] = func(x);%,[]);

% Initialize search directions
r = -g;
s = precon_func(r);
d = s;

delta_new = dot(r,d);
delta0 = delta_new;

stimfilt_idx = 2:1+basepars.n*basepars.Mk;

% initial filtering of stimulus
%kx_base = filterstimulus_train2_nonsep(x,basepars,stimpars,trainpars);

1;
npars = 1+basepars.n*basepars.Mk + basepars.nofilters_postspike + (basepars.Nneurons-1)*basepars.nofilters_coupling;

while (i < pars.max_iter && delta_new > pars.tol^2*delta0)

   1; 
   j = 0;
   fprintf('CG iteration %d: f=%f\t g=%f\n',i,f,norm(g));

   % prefilter the stimulus along the search direction
   %fprintf('done.\nFiltering in search direction...');
   xtest = x;   
   xtest(stimfilt_idx) = d(stimfilt_idx);
   %kx_dir = filterstimulus_train2_nonsep(xtest,basepars,stimpars,trainpars);
   %fprintf('done.\n');
   
   delta_d = dot(d,d);
   
   % Do the linesearch along direction d with specified method
   alpha_total = 0; % cumulative amount moved in direction d during linesearch
   switch pars.ls_method
    
       case 'secant' % Secant linesearch
           alpha = - pars.step;
           [fstep gstep] = func(x + pars.step*d);%,kx_base+pars.step*kx_dir);
           eta_prev = dot(gstep,d);
           alpha_total = 0; 
           while 1
               if (mod(j,10) == 0)
               % fprintf('\t Secant Linesearch iteration %d\n',j);
               end
               eta = dot(g,d);
               alpha = alpha*eta/(eta_prev-eta);
               x = x + alpha*d;
               alpha_total = alpha_total + alpha;
               [f g] = func(x);%,kx_base + alpha_total*kx_dir);
               eta_prev = eta;
               j = j+1;
               
               if (j >= pars.max_ls_iter || alpha^2*delta_d <= pars.ls_tol^2)
                   break;
               end
           end
       case 'newton' % Newton-Raphson linesearch
           1;
           while 1
               1;
               if (mod(j,10) == 0)
               % fprintf('\t Newton Linesearch iteration %d\n',j);
               end
               [f g cifs] = func(x);%,kx_base + alpha_total*kx_dir);
               1;
               dHd = 0;
               if (basepars.Nneurons > 1)
                   fprintf('ERROR: newton linesearch not configured for multiple neurons yet!\n');
               end
               for k=1:basepars.Nneurons
                idx = (k-1)*npars+1:k*npars;
                dHd = dHd - dot(d(idx),multres_product([ones(1,size(stimpars.x,2));trainpars.lgrad{k}],trainpars.psbasisGrad{k},-trainpars.dt.*cifs(:,k),d(idx(1:1+basepars.n*basepars.Mk)),d(idx(1+basepars.n*basepars.Mk+1:end))));
               end
               alpha = - dot(g,d)/dHd;
               x = x + alpha*d;
               alpha_total = alpha_total + alpha;
               j = j+1;
               if (j >= pars.max_ls_iter || alpha^2*delta_d <= pars.ls_tol^2)
                   break;
               end
           end
   end
   
   if (j == pars.max_ls_iter)
       fprintf('\tLinesearch: Max ls iters achieved!\n');
   else
       fprintf('\t Linesearch: achieved appropriate step after %d iterations. Step norm is %f\n',j,alpha_total.*norm(d));
   end
   
   r = -g;
   delta_old = delta_new;
   delta_mid = dot(r,s);
   s = precon_func(r);
   delta_new = dot(r,s);

   %fprintf('Val/Tol is %f/%f\n',delta_new,pars.tol^2*delta0);
   beta = (delta_new-delta_mid)/delta_old;

   % Update the prefiltered stimulus
   %kx_base = kx_base + kx_dir*alpha_total;
   
   k = k+1;
   
   if (k == length(x) || beta <= 0) % Reset the search if necessary
       d = s;
       k = 0;
   else
       d = s + beta*d;
   end
   
   i = i+1;
   
   if (isfield(pars,'plotfcn')) % Plot any suppled functions at every iteration
       pars.plotfcn(x);
   end
   
end 
[f g] = func(x);%,[]);
fprintf('Converged to solution after %d iterations. Gradient norm is %f!\n',i,norm(g));

xstar = x;