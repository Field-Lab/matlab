function [ll_grad,ll_hess] = train_ll_grad4_lite(p,basepars,~,trainpars,lcifs,cifs)
% Function that computes the GRADIENT of the log likelihood of the data D for the parameters p

%fprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
%keyboard

hessFlag = 0; %( ~strcmp(basepars.hessmode,'mult') || ~strcmp(basepars.filtermode,'sta') );
% (~(1) || ~(1)) = (0 || 0) = 0

fprime = basepars.Nprime(lcifs);
%kernel = fprime./cifs.*full(double(trainpars.D(:,trainpars.baseneuron_idx))) - trainpars.dt.*fprime; % approx. for dt small
kernel = fprime./cifs.*full(double(trainpars.D(:))) - trainpars.dt.*fprime; % approx. for dt small

ll_grad = multisample_gradHess_structs_lite(basepars,1,kernel,[],trainpars);

if (nargout > 1)
   fdprime = basepars.Ndoubleprime(lcifs);
   %kernel2 = (fdprime.*cifs - fprime.^2)./(cifs.^2).*full(double(trainpars.D(:,trainpars.baseneuron_idx))) - fdprime.*trainpars.dt;
   kernel2 = (fdprime.*cifs - fprime.^2)./(cifs.^2).*full(double(trainpars.D(:))) - fdprime.*trainpars.dt;
   %keyboard
   if (hessFlag)
      ll_hess = zeros(length(p),length(p));
      %keyboard
   end
else
   kernel2 = [];
end

% Make a numpars x pars.maxt matrix where the jth column is the
% gradient of the input term at time j for this trial (for the nth neuron)


%keyboard
if (nargout >1 && ~hessFlag)
   %    fprintf('Creating Hinfo structure\n');
   ll_hess = struct('kernel',kernel,'kernel2',kernel2);
end
