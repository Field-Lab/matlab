% Training the parameters of the glm to data given a fixed set of stimulus filters

% Parameters to train:
% 1. stimulus filter (space x time) for each neuron
% 2. the coupling filters h (for each pair of neurons) - this is represented as cooridinates wrt some finite basis (of exponentials).
%    This is a 3D matrix (basis-coordinate-vector x neuron 1 x neuron 2)
% 3. the base firing rate of each neuron. This is a vector (neurons x 1)

% The function to be optimized is the log-likelihood of the data:

% sum_i[ sum_{t_{i,j}}[ log(L_{i,t})] - sum_t[L_{i,t}dt] ]


% the arguments:
% p0 - initial search point - the vector is blocked according to neurons.
%    - within each neuron block it is ordered as:
%    - base firing rate, stimulus filter, postspike basis coeffs, coupling basis coeffs (ordered by neuron)

% basepars - structure with the following fields:
    % maxt - # of timebins
    % n - spatial dimension
    % Nneurons - number of neurons 
    % postspike_basis - matrix where each column is a postspike basis function
    % couping_basis - "" coupling basis fns
    % nofilters_postspike - number of postspike basis fns
    % nofilers_coupling - "" coupling
    % fac - ratio of stimulus temporal resolution to spike train resolution
    %       - IMPORTANT!
    % fratio - ratio of timespan of postspike filter to stimulus filter
    % 
    
% stimpars - structure with the following fields:
    % x - stimulus (space x time)
    % dt - temporal resolution of stimulus
    
% trainpars - structure with the following fields:
    % dt - temporal resolution of spike train time bins
    % D - sparse logical matrix where ith column is the spike train of the ith neuron
    
% gopts - optimization options to pass to fminunc



% Returns:
% pstar - optimal solution returned by fminunc
% fstar, gstar,Hstar - function, gradient, and hessian value at the optimal point
% eflag,output - exit flag and output structure returned by fminunc

function [pstar fstar gstar Hstar eflag output] = train_glm2(p0,basepars,stimpars,trainpars,gopts,NITER)
% fprintf('check point @ train_glm2 (edoi)')
% keyboard

if (~exist('NITER','var'))
    NITER = 1;
end
fl_filt_mode = ( strcmp(basepars.filtermode,'nonsep') || strcmp(basepars.filtermode,'raw'));

% Initialize the trainpars structure

% It contains any values that need not be computed on every optimization
% iteration such as the hessian spacetime constant matrix (i.e., the spike trains, temporally
% filtered stimulus), the convolutions of the spike trains with the
% postspike/coupling basis functions
1;
%if (~isempty(basepars.ktime_basis) && ~isfield(trainpars,'hst'))
%    trainpars.hst = hess_spacetimeconst(stimpars.x,basepars.ktime_basis,stimpars.dt);
%end

%if (isempty(basepars.ktime_basis) && ~isfield(trainpars,'hstraw'))
%    trainpars.hstraw = hess_spacetimeconst_raw(stimpars.x,basepars.Mk);
%end


% There is a slight shift in the alignment of spike train recording and
    % stimulus presentation  - Jon Shlens said it is approx. half of the f
    % ratio +/- 1 in the fine representation
%basepars.ps_shift_val = ceil(f/2 + 1);

if (~isfield(trainpars,'psbasisGrad'))
    fprintf('Computing ps basis grad...');
    trainpars.psbasisGrad = grad_basis(trainpars.D(:,1:basepars.Nneurons),basepars.postspike_basis); 
    fprintf('done.\n');
    %trainpars.psbasisGrad = trainpars.psbasisGrad(:,ps_
end

if (size(trainpars.D,2) > 1 && ~isfield(trainpars,'cpbasisGrad'))
    fprintf('Computing cp basis grad...');
    trainpars.cpbasisGrad = grad_basis(trainpars.D,basepars.coupling_basis);
    fprintf('done.\n');
elseif (size(trainpars.D) == 1)
    trainpars.cpbasisGrad = [];
end


if ( fl_filt_mode && ~isfield(trainpars,'lgrad') )
   fprintf('Computing lgrad for nonsep mode...');
   trainpars.lgrad = nonsep_lgrad(basepars,stimpars,trainpars); % edoi
   fprintf('done.\n');
end


% Optimization use fminunc

if (~exist('gopts','var') || isempty(gopts)) % Set to default options
    fprintf('USING DEFAULT TRAINING OPTIONS!\n');
    gopts = optimset('GradObj','on','Hessian','on','TolFun',10^(-6),'MaxIter',200,'display','iter');%,'PlotFcns',@optimplotfval);%,'DerivativeCheck','on')
end


%-- plot params during iteration
% if (isempty(optimget(gopts,'PlotFcns'))) % plotting settings if not already set
%     pfcn = setup_plotfcns(basepars.filtermode,basepars,stimpars,trainpars);
%     gopts = optimset(gopts,'PlotFcns',pfcn);
% end

1; % Here is where you want to validate your gradient and hessian functions
% disable plotting
if 0
%vartest = [1 56 8000 11525 11533];
%vartest = [1 100 500 1200 1450 2000 2050 2500 2700 2890];
%check_gradhess(0.1.*randn(size(p0)),@(p) ll_func2(p,basepars,stimpars,trainpars), @(Hinfo,Y) ll_hess_mult(Hinfo,Y,basepars,stimpars,trainpars), vartest);
vartest = [1 5 10 20 25 35 50];
check_gradhess(p0,@(p) ll_func2(p,basepars,stimpars,trainpars), [], vartest);
end
%gopts = optimset(gopts,'PlotFcns',[]);

if ( (strcmp(basepars.hessmode ,'mult') && fl_filt_mode && (~isfield(basepars,'normconstr') || basepars.normconstr == 0)) ) || ...
   ( (strcmp(basepars.hessmode ,'mult') && strcmp(basepars.filtermode,'raw') && (~isfield(basepars,'normconstr') || basepars.normconstr == 0)) )
    gopts = optimset(gopts,'HessMult',@(Hinfo,Y) ll_hess_mult(Hinfo,Y,basepars,stimpars,trainpars));
    pstar = p0;
    1;
%    fprintf('check-2')
%    keyboard
    for j=1:NITER
        
        if (isfield(basepars,'ext_timepts') && ~isempty(basepars.ext_timepts))
            mtx = sum(interpmtx(basepars.interpfilt,basepars.ext_timepts,basepars.maxt));
            A = zeros(1,length(pstar));
            A(get_pars_idx(basepars,1,size(trainpars.D,2),'ext')) = mtx;
            fprintf('Constraining extrinsic signal to have mean 0!\n');           
            [pstar fstar eflag output] = fmincon(@(p)ll_func2(p,basepars,stimpars,trainpars),pstar,[],[],A,0,[],[],[],gopts); 
        else
            [pstar fstar eflag output] = fminunc(@(p)ll_func2(p,basepars,stimpars,trainpars),pstar,gopts);
        end
    end
      
    [fstar gstar] = ll_func2(pstar,basepars,stimpars,trainpars);
    [Hstar] = [];%ll_hess_mult(Hinfo,eye(size(gstar,1)),basepars,stimpars,trainpars);
    
elseif ( (strcmp(basepars.hessmode ,'mult') && fl_filt_mode && (isfield(basepars,'normconstr') && basepars.normconstr)) ) || ...
   ( (strcmp(basepars.hessmode ,'mult') && strcmp(basepars.filtermode,'raw') && (isfield(basepars,'normconstr') && basepars.normconstr)) )
    % Stimulus filter norm-constraints in nonseparable mode 
    gopts_constr = optimset(gopts,'Hessian','user-supplied','HessMult',@(x,lambda,Y) ll_hess_mult_normconst(x,lambda,Y,basepars,stimpars,trainpars),'GradConstr','on','Algorithm','interior-point','SubproblemAlgorithm','cg');
    [pstar fstar] = fmincon(@(p)ll_func2(p,basepars,stimpars,trainpars),p0,[],[],[],[],[],[],@(x) norm_constr(x,basepars,stimpars,trainpars),gopts_constr);    
    [fstar gstar] = ll_func2(pstar,basepars,stimpars,trainpars);
    Hstar = [];
elseif (strcmp(basepars.filtermode,'sep_raw') && (isfield(basepars,'normconstr') && basepars.normconstr))
    % Spatial filter norm constraints in separable mode    
    fprintf('CONSTRAINING NORMS OF SPATIAL FILTERS....\n');
    gopts_constr = optimset(gopts,'Hessian','user-supplied','HessFcn',@(x,lambda) hessian_sep_normconstr(x,lambda,basepars,stimpars,trainpars),'GradConstr','on','Algorithm','interior-point','SubproblemAlgorithm','cg');
    [pstar fstar] = fmincon(@(p)ll_func2(p,basepars,stimpars,trainpars),p02,[],[],[],[],[],[],@(x) norm_constr(x,basepars,stimpars,trainpars),gopts_constr);    
    [fstar gstar Hstar] = ll_func2(pstar,basepars,stimpars,trainpars);

else
    outflag = 0;
    counter = 0;
    while (~outflag)
        pstar = p0;
        1;
        
        %check_gradhess(0.1.*randn(size(p0)),@(p) ll_func2(p,basepars,stimpars,trainpars), []);
        
        
        % Put a norm 1 constraint on the spatial filters
        if 0
        basepars.normconstr_p = 1;
        basepars.normconstr = 1;
        nc = @(x) norm_constr(x,basepars,stimpars,trainpars);
        gopts_constr = optimset(gopts,'Hessian','user-supplied','HessFcn',@(x,lambda) hessian_sep_normconstr(x,lambda,basepars,stimpars,trainpars),'GradConstr','on','Algorithm','interior-point','SubproblemAlgorithm','cg');
        pstar = normalize_sep_filts(pstar,basepars,size(trainpars.D,2),basepars.normconstr_p);
        fprintf('Constraining Lp norm of spatial filters to be less than 1!\n');
 
        end
        
        
        % Linear bounds on the spatial filers
        if 0
        %lb = -inf.*ones(size(pstar));
        %ub = inf.*ones(size(pstar));
        %spaceidx = get_pars_idx(basepars,1,size(trainpars.D,2),'kspace');
        %lb(spaceidx) = -1;
        %ub(spaceidx) = 1;
        %pstar = bound_sep_filts(pstar,basepars,size(trainpars.D,2));
        end
        
        % Put a norm constraint on the PS filters
        if (isfield(basepars,'psnormconstr') && basepars.psnormconstr > 0)
            fprintf('Constrining the norm of the PS filter coefficients.\n');
            nc = @(x) ps_norm_constr(x,basepars,stimpars,trainpars);
            gopts_constr = optimset(gopts,'Hessian','user-supplied','HessFcn',@(x,lambda) hessian_sep_psnormconstr(x,lambda,basepars,stimpars,trainpars),'GradConstr','on','Algorithm','interior-point','SubproblemAlgorithm','cg');
        end
        
        %fprintf('Adjusting extrinsic signal if necessary.\n');
        pstar = adjust_baserate_ext(pstar,basepars,stimpars,trainpars);
        if (counter > 0)
            pstar = 0.1.*randn(size(pstar)); 
            gopts = optimset(gopts,'MaxIter',120); % allow more iterations
        end
        
        1;
        
        
        for j=1:NITER
            
            if (isfield(basepars,'ext_timepts') && ~isempty(basepars.ext_timepts))
                
                mtx = sum(interpmtx(basepars.interpfilt,basepars.ext_timepts,basepars.maxt));
                A = zeros(1,length(pstar));
                A(get_pars_idx(basepars,1,size(trainpars.D,2),'ext')) = mtx;
                fprintf('Constraining extrinsic signal to have mean 0!\n');
                if (~exist('nc','var'))%
                    [pstar fstar eflag output] = fmincon(@(p) ll_func2(p,basepars,stimpars,trainpars), pstar, [], [], A, 0, [], [], [], gopts);
                    % [pstar fstar eflag output] = fmincon(@(p) ll_func2(p,basepars,stimpars,trainpars), pstar, [], [], A, 0, lb, ub, [], gopts);
                else
                    [pstar fstar eflag output] = fmincon(@(p) ll_func2(p,basepars,stimpars,trainpars), pstar, [], [], A, 0, [], [], nc, gopts_constr);
                end
                
            else
                if (~exist('nc','var'))
                    %[pstar fstar eflag output] = fmincon(@(p) ll_func2(p,basepars,stimpars,trainpars), pstar, [], [], [], [], lb, ub, [], gopts);
                    fprintf('unconstrained optimization..')
                    [pstar fstar eflag output] = fminunc(@(p) ll_func2(p,basepars,stimpars,trainpars),pstar,gopts);
                else
                    fprintf('constrained optimization..')
                    [pstar fstar eflag output] = fmincon(@(p) ll_func2(p,basepars,stimpars,trainpars), pstar, [], [], [], [], [], [], nc, gopts_constr);
                end
                
            end
        end
        
        if 0%(eflag == -1)
            outflag = 0;
            counter = counter + 1;
            fprintf('Optimization terminated incorrectly. restarting...\n');
        else
            outflag = 1;            
        end
        
        1;
        if 0
            % Fix the shape of the filter and optimize
            K = get_model_filters(pstar,1,basepars,1);
            pstarfixfilt = [pstar(get_pars_idx(basepars,1,size(trainpars.D,2),'ext'));
                pstar(get_pars_idx(basepars,1,size(trainpars.D,2),'b'))
                1;
                pstar(get_pars_idx(basepars,1,size(trainpars.D,2),'ps'))];
            basepars2 = basepars; basepars2.filtermode = 'fixfilt';
            basepars2.K{1} = K;
            for j=1:basepars.Nneurons
                trainpars.kx = stimpars.dt .* sum(fastconv(stimpars.x(basepars.crop_idx(:,j),:),basepars2.K{j},length(basepars.crop_idx(:,j)),size(stimpars.x,2)))';
            end
            pfcn2 = setup_plotfcns(basepars2.filtermode,basepars2,stimpars,trainpars);
            gopts2 = optimset(gopts,'PlotFcns',pfcn2);
            [pstarfixfilt] = fminunc(@(p) ll_func2_fixfilt(pstarfixfilt,basepars2,stimpars,trainpars), pstarfixfilt, gopts2);
        end
        
        
        [fstar gstar Hstar] = ll_func2(pstar,basepars,stimpars,trainpars);
    end
end
