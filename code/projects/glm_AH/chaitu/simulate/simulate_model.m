% Function to simulate the GLM using paramters specified in vector p (and
% basepars) driven by a stimulus X for ntrials

% Returns the binary matrix of spikes and the log of the resulting CIFs (time x neuron x trial)

function [spikes lcifs kx] = simulate_model(basepars,p,X,Deff,dt,ntrials,opt_pars)

if (~isempty(X))
   basepars.maxt = size(X,2);
elseif (exist('opt_pars','var'))
   basepars.maxt = floor(1.1*3/2*floor(opt_pars.period/(dt*basepars.fac)));
end


% Number of time bins in spike train/cif

spT = basepars.maxt * basepars.fac;

% Initalize output
lcifs = zeros(spT,basepars.Nneurons,ntrials);
spikes = (logical(false(spT,basepars.Nneurons,ntrials)));

% Setup a stimpars struct just to filter it
stimpars.x = X;
stimpars.dt = dt*basepars.fac;

% Filter the stimulus
if (~exist('opt_pars','var'))
   % Filter the stimulus
   if (isfield(basepars,'XsqK') && basepars.XsqK)
      [kx kxsq] = filterstimulus_train(p,basepars,stimpars,struct('D',[0]));
      kx = kx + kxsq;
   else
      %keyboard
      %if strcmp(basepars.filtermode,'lcs1')
      %   test_sp_gen_sig_lcs1;
      %end
      kx = filterstimulus_train(p,basepars,stimpars,struct('D',[0]));
   end
end

N = basepars.Nneurons;
Neff = N + size(Deff,2);

numpars_perneuron = get_npars(basepars,Neff);
offset = get_pars_idx(basepars,1,Neff,'ps');
offset = offset(1)-1;

% Recover indices for the postspike filters
add_term = repmat(numpars_perneuron,basepars.nofilters_postspike,N)*diag(0:(N-1)); % basepars.nofilters_postspike x N
idx = repmat((offset+1:offset+basepars.nofilters_postspike)',1,N) + add_term; % basepars.nofilters_postspike x N
ps_weights = reshape(p(idx(:)),basepars.nofilters_postspike,N);
%PS = basepars.postspike_basis*ps_weights; % basepars.Mhist x N
%[mhist nps] = size(PS);
offset = offset + basepars.nofilters_postspike;

if (Neff > 1)
   % Coupling filters corresponding to these parameters (dim: basepars.Mhist x neuron1*neuron2)
   add_term = repmat(numpars_perneuron,(Neff-1)*basepars.nofilters_coupling,N)*diag(0:N-1); % (N-1)*nofitlers_coupling x N
   idx = repmat((offset+1:offset+(Neff-1)*basepars.nofilters_coupling)',1,N) + add_term; % (N-1)*basepars.nofilters_coupling x N
   cp_weights = reshape(p(idx(:)),basepars.nofilters_coupling,N*(Neff-1)); % basepars.nofilters_coupling x N*(N-1)
   %    CP = basepars.coupling_basis*cp_weights;
   %    ncoupling = size(c_weights,1);
else
   cp_weights = [];
end

1;




exflag = 1;
1;
MAX_ITER = 10;

%disppercent(-inf,'Starting simulation.');
for j=1:ntrials
   
   % set a seed for exprnd
   %rand('twister',j);
   
   counter = 1;
   fprintf('Trial %d/%d\n',j,ntrials);
   %if (mod(j,10) == 0)
   %disppercent(j/ntrials);
   %fprintf('Completed %d trials.\n',j);
   %end
   
   while(exflag && counter < MAX_ITER)
      
      if (counter > 1)
         fprintf('Restarting trial %d: iteration %d\n',j,counter);
      end
      
      1;
      if (exist('opt_pars','var'))
         %Regenerate random stimulus
         switch (opt_pars.mode)
            
            case 'gwn'
               stimpars.x = opt_pars.mu + opt_pars.sig.*randn(basepars.n,basepars.maxt);
               
            case 'switch'
               
               nframes_perswitch = floor((opt_pars.period/2)/(dt*basepars.fac));
               
               switch(opt_pars.noisetype)
                  case 'bwn'
                     
                     
                     
                     idx = 1:nframes_perswitch;
                     if (isfield(stimpars,'x'))
                        stimpars = rmfield(stimpars,'x');
                     end
                     
                     if (strcmp(basepars.stim_type,'lv'))
                        1;
                        %dim = basepars.n/(4^2);
                        %stimpars.x(:,idx) = opt_pars.mu2-opt_pars.sig2+2.*opt_pars.sig2.*(rand(dim,length(idx))<0.5);
                        %idx = nframes_perswitch+1:2*nframes_perswitch;
                        %stimpars.x(:,idx) = opt_pars.mu1-opt_pars.sig1+2.*opt_pars.sig1.*(rand(dim,length(idx))<0.5);
                        %idx = 2*nframes_perswitch+1:basepars.maxt;
                        %stimpars.x(:,idx) = opt_pars.mu2-opt_pars.sig2+2.*opt_pars.sig2.*(rand(dim,length(idx))<0.5);
                        %stimpars.x = resize_stimulus(stimpars.x,[sqrt(dim) sqrt(dim)],[sqrt(basepars.n) sqrt(basepars.n)]);
                        %[s1 t1 s2 t2] = get_sep_filters(p,basepars.nofilters_kspace,basepars.Mk,basepars.kspace_basis,eye(basepars.Mk),0);
                        1;
                        kx = create_random_upsampled_kx(basepars,p,basepars.maxt,stimpars.dt,opt_pars);
                     else
                        dim = basepars.n;
                        stimpars.x(:,idx) = opt_pars.mu2-opt_pars.sig2+2.*opt_pars.sig2.*(rand(dim,length(idx))<0.5);
                        idx = nframes_perswitch+1:2*nframes_perswitch;
                        stimpars.x(:,idx) = opt_pars.mu1-opt_pars.sig1+2.*opt_pars.sig1.*(rand(dim,length(idx))<0.5);
                        idx = 2*nframes_perswitch+1:basepars.maxt;
                        stimpars.x(:,idx) = opt_pars.mu2-opt_pars.sig2+2.*opt_pars.sig2.*(rand(dim,length(idx))<0.5);
                        kx = filterstimulus_train(p,basepars,stimpars,struct('D',[0]));
                     end
                     
               end
               
            case 'switch-supplied'
               stimpars.x = X(:,:,j);
               
         end
      end
      
      bidx = get_pars_idx(basepars,1,Neff,'b');
      bvals = p(bidx(1):numpars_perneuron:length(p))';
      1;
      if (isempty(Deff))
         [exflag sp lcifs(:,:,j)] = simulate_trains_fast4(basepars,dt,kx,bvals,[],ps_weights,cp_weights);
      else
         [exflag sp lcifs(:,:,j)] = simulate_trains_fast4(basepars,dt,kx,bvals,Deff(:,:,j),ps_weights,cp_weights);
      end
      counter = counter+1;
   end
   
   if (exflag)
      spikes = [];
      fprintf('Failed to simulate with this model.\n');
      return;
   end
   exflag = 1;
   
   [I J] = find(sp);
   for k=1:length(I)
      spikes(I(k),J(k),j) = true;
   end
end
%disppercent(inf);