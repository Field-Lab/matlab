% It returns function, gradient (and Hessian) value of the log likelihood of the data.
% lite: edoi, 2012-01-10.

% t - time vector (unit of time: )
% maxt, n, Nneurons - time, spatial, neuronal dimensions
% x - stimulus
% D - data (3D binary matrix: time x trial x neuron)
% H is a dummy variable in lite version.

function [f g H lcifs] = ll_func2_lite(p,basepars,stimpars,trainpars,kx_opt)

%hessFlag = ( ~strcmp(basepars.hessmode,'mult') || ~strcmp(basepars.filtermode,'sta') );
hessFlag = 0;

%-- If any parameters are frozen at designated values, reset them to their designated values
frozenflag = isfield(basepars,'frozen_idx') && isfield(basepars,'frozen_vals');
if (frozenflag)
   p(basepars.frozen_idx) = basepars.frozen_vals;
end

%-- calculate the log likelihood of the data
if ( strcmp(basepars.filtermode,'raw') || strcmp(basepars.filtermode,'rk2') )
   %fprintf('\nassuming the model is GLM\n')
   [logprob lcifs cifs] = train_ll3_glm(p,basepars,stimpars,trainpars);
elseif (exist('kx_opt','var') && ~isempty(kx_opt))
   %fprintf('\nassuming the model is vGLM\n')
   [logprob lcifs cifs] = train_ll3_lite(p,basepars,stimpars,trainpars,kx_opt);
else
   %fprintf('\nassuming the model is vGLM\n')
   [logprob lcifs cifs] = train_ll3_lite(p,basepars,stimpars,trainpars);
end
f = full(-logprob); % we want to minimize negative log-likelihood


%-- gradient
if nargout == 1
   1;
elseif nargout < 3
   g = train_ll_grad4_lite(p,basepars,stimpars,trainpars,lcifs,cifs);
else
   [g H] = train_ll_grad4_lite(p,basepars,stimpars,trainpars,lcifs,cifs);
   if (hessFlag)
      H = -H;
   fprintf('####')
   end
   %keyboard
   %fprintf('Condition number of total hessian is %f\n',cond(H));
end


if exist('g','var')
   if (frozenflag)
      g(basepars.frozen_idx) = 0;
   end
   g = -g;
   
%   fprintf('f=%f    gnorm=%f\n',f,norm(g));
end
