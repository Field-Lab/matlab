function [xstar fstar gstar] = train_glm_alternate(x0,basepars,stimpars,trainpars,gopts)



% Preconditioned nonlinear CG

% First lets setup the preconditioner function - we need to first calculate
% the inverse of n*E[xx'exp(b + x'*theta)]= n*exp(b + theta'*S*theta/2 + mu'*theta).*(S + (mu+S*theta)*(mu + S*theta)');
% where S = E[xx'], mu = E[x]

1;

if 0
switch(basepars.noise_type)
    case'gwn'
        % Gaussian white noise case
        %fprintf('Computing the covariance and its inverse...');
        %mu = stimpars.dt * repmat(0.5,basepars.n*basepars.Mk,1);
        %sig = stimpars.dt * sqrt(1/2 * ( (0.02-0.5)^2 + (0.98-0.5)^2));
        %S = sig^2 * eye(basepars.n*basepars.Mk) + mu*mu'; % Init the covariance matrix
        %O = fft(eye(basepars.n*basepars.Mk)); % fourier basis matrix
        %d = fft(S(1,:))'.*size(S,2);
        %row1 = ifft(1./d);
        %Sinv = size(S,2) .* toeplitz([row1(1);row1(end:-1:2)],row1);
        %fprintf('done.\n');
        %Sinv = mult_col(real(O),1./d)*real(O)' + mult_col(imag(O),1./d)*imag(O)';
        
    case 'bwn'
        % Binary white noise case
        a = stimpars.dt * (0.98-0.02);
        b = stimpars.dt * 0.02;
    case 'g1fn'
        % Gaussian 1/f noise case
        fprintf('Computing the covariance and its inverse...');
        mu_val = 0.5;
        sig_val = sqrt(0.025);        
        [S Sinv] = compute_1f_covariance2(40,2,320,160,5,basepars.n,basepars.crop_idx,12,basepars.Mk,mu_val,sig_val);
        mu = mu_val.*ones(basepars.n*basepars.Mk,1);
        % Scale the mean and covariance appropriately
        S = (stimpars.dt^2) .* S;
        Sinv = Sinv ./ (stimpars.dt^2);
        mu = stimpars.dt .* mu;
end
end

% Initialize
xstar = x0;

% Set up the CG pars
pars.max_outer_iter = 5;
pars.max_iter = 10;
pars.tol = 10^(-3);
pars.step = 10^(-8);
pars.max_ls_iter = 100;
pars.ls_tol = 10^(-8);
pars.plotfcn = @(x) cg_plotpars(x,basepars,stimpars,trainpars);
pars.ls_method = 'secant';


1;

func = @(x) ll_func_cg(x,basepars,stimpars,trainpars);
%func = @(x) ll_func_cg(x,basepars,stimpars,trainpars);

1;

fprintf('Starting nonlinear CG...\n');

for j=1:pars.max_outer_iter
    
    fprintf('\n\nCG Outer iteration number %d...\n',j);
    % Use the Woodbury matrix lemma to estimate the Hessian inverse
    theta = xstar(2:1+basepars.n*basepars.Mk); % mn x 1
 
    fprintf('Computing the preconditioner...'); % Use woodbury lemma
    
    Hinv_partial_allcells = zeros(basepars.n*basepars.Mk,basepars.n*basepars.Mk,basepars.Nneurons);
    
    for k=1:basepars.Nneurons
        
        if 0
            switch (basepars.noise_type)
            
            case 'gwn'
                
                % Gaussian white noise case
                %mst = mu + S*theta;
                %c_const = max(10^(-6),exp(-xstar(1) - theta'*S*theta/2 - mu'*theta));
                %Hinv_partial =  c_const .* (Sinv - ((Sinv*(mst))*((mst)'*Sinv))./(1+ (mst'*(Sinv*mst))));
                
            case 'bwn'
                1;
                % Binary white noise case
                d = length(theta);
                eatheta = exp(a.*theta);
                ctheta = -trainpars.dt .* exp(xstar(1) +d*log(0.5) + b*sum(theta) + sum(log(1+eatheta))); % compute in logspace to avoid numerical errors
                ctheta = exp(xstar(1)) * (0.5)^(d) * exp(b*sum(theta)) * prod(ones(d,1) + eatheta);
                alpha = eatheta./(1+eatheta);
                abvec = (a.*alpha + b.*ones(d,1));
                M1diag = a^2.*(alpha - alpha.^2);
                H_partial = - size(stimpars.x,2)*basepars.fac*ctheta*(diag(M1diag) + abvec*abvec');
                M1inv = diag(1./M1diag);
                Hinv_partial_allcells(:,:,k) = -1/(ctheta*size(stimpars.x,2)*basepars.fac) .* (M1inv - (M1inv*abvec)*(abvec'*M1inv)./(1+dot(abvec,M1inv*abvec)));
                
                
            case 'g1fn'
                %mst = mu + S*theta;
                %Hinv_partial1 = (Sinv - ((Sinv*(mst))*((mst)'*Sinv))./(1+ (mst'*(Sinv*mst))));
                %ctheta = 1/(trainpars.dt .* size(stimpars.x,2)*basepars.fac) .* exp(-xstar(1) - theta'*S*theta/2 - mu'*theta);
                %Hinv_partial = Hinv_partial1.* ctheta;
                Hinv_partial_allcells(:,:,k) =  exp(-xstar(1)) * 1/(size(stimpars.x,2)*basepars.fac) * Sinv;
                % Others...
                
            end
        end
    end
    
    % Compute the inverse
    %fprintf('Computing the hessian of the other terms...');
    %tic;
    %[f g cifs] = func(xstar,[]);
    %PSH = trainpars.dt .* mult_col(trainpars.psbasisGrad,cifs)*trainpars.psbasisGrad';
    %fprintf('done.\n');
    %toc;
    %Hinv = blkdiag(eye(1),Hinv_partial,eye(basepars.nofilters_postspike));
    %Hinv = blkdiag(1/(trainpars.dt.*sum(cifs)),Hinv_partial,inv(PSH));
    precon_func = @(x) x;%blkdiag_multiply(Hinv_partial_allcells,x,2:1+basepars.n*basepars.Mk); %Hinv*x;
    
    fprintf('done.\n');
   
    tic; xstar = pcnlcg(xstar,func,pars,basepars,stimpars,trainpars,precon_func); toc;
    
end
1;
[fstar gstar] = func(xstar);%,[]);


%% 
if 0 
    gopts2 = optimset(gopts,'HessMult',@(Hinfo,Y) ll_hess_mult(Hinfo,Y,basepars,stimpars,trainpars),'PrecondBandWidth',100);
    pfcn = setup_plotfcns(basepars.filtermode,basepars,stimpars,trainpars);
    gopts2 = optimset(gopts2,'PlotFcns',pfcn)
    gopts2 = optimset(gopts2,'MaxIter',2)
    tic; [pstar2 f2 g2 h2] = train_glm2(xstar,basepars,stimpars,trainpars,gopts2); toc;
end