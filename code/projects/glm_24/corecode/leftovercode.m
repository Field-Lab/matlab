%%%%%%% Train_glm2_AH   un used calls for fullr ank raw or seperable filters%%%%%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%if (~isempty(Basepars.ktime_basis) && ~isfield(Trainpars,'hst'))
%    Trainpars.hst = hess_spacetimeconst(Stimpars.movie_ROI,Basepars.ktime_basis,Stimpars.dt);
%end

%if (isempty(Basepars.ktime_basis) && ~isfield(Trainpars,'hstraw'))
%    Trainpars.hstraw = hess_spacetimeconst_raw(Stimpars.movie_ROI,Basepars.Mk);
%end


% There is a slight shift in the alignment of spike train recording and
    % stimulus presentation  - Jon Shlens said it is approx. half of the f
    % ratio +/- 1 in the fine representation
%Basepars.ps_shift_val = ceil(f/2 + 1);


%%% AH not really sure what lgrd means here
if ( filtermode && ~isfield(Trainpars,'lgrad') )
   fprintf('Computing lgrad for nonsep mode...');
   Trainpars.lgrad = nonsep_lgrad(Basepars,Stimpars,Trainpars);  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WORK HERE %%%%%%%%%%%%%%%%%%%%%%%%%%
   fprintf('done.\n');
end

%%%% AH checking  for params gopt for the fminunc
if (~exist('gopts','var') || isempty(gopts)) % Set to default options
    fprintf('USING DEFAULT TRAINING OPTIONS!\n');
    gopts = optimset('GradObj','on','Hessian','on','TolFun',10^(-8),'MaxIter',50,'display','iter');%,'PlotFcns',@optimplotfval);%,'DerivativeCheck','on')
end


%-- plot params during iteration
% if (isempty(optimget(gopts,'PlotFcns'))) % plotting settings if not already set
%     pfcn = setup_plotfcns(Basepars.filtermode,Basepars,Stimpars,Trainpars);
%     gopts = optimset(gopts,'PlotFcns',pfcn);
% end

% Here is where you want to validate your gradient and hessian functions
% disable plotting
if 0
%vartest = [1 56 8000 11525 11533];
%vartest = [1 100 500 1200 1450 2000 2050 2500 2700 2890];
%check_gradhess(0.1.*randn(size(p0)),@(p) ll_func2(p,Basepars,Stimpars,Trainpars), @(Hinfo,Y) ll_hess_mult(Hinfo,Y,Basepars,Stimpars,Trainpars), vartest);
vartest = [1 5 10 20 25 35 50];
check_gradhess(p0,@(p) ll_func2(p,Basepars,Stimpars,Trainpars), [], vartest);
end
%gopts = optimset(gopts,'PlotFcns',[]);



%%% mult means something about being quasi newton.. not quite sure


%%% IF Raw filter and no norm contraint has yet to be put in place
if ( (strcmp(Basepars.hessmode ,'mult') && filtermode && (~isfield(Basepars,'normconstr') || Basepars.normconstr == 0)) )  
    gopts = optimset(gopts,'HessMult',@(Hinfo,Y) ll_hess_mult(Hinfo,Y,Basepars,Stimpars,Trainpars));
    pstar = p0;
    1;
    % fprintf('check-2');    keyboard
    for j=1:NITER
        if (isfield(Basepars,'ext_timepts') && ~isempty(Basepars.ext_timepts))
            mtx = sum(interpmtx(Basepars.interpfilt,Basepars.ext_timepts,Basepars.maxt));
            A = zeros(1,length(pstar));
            A(get_pars_idx(Basepars,1,size(Trainpars.D,2),'ext')) = mtx;
            fprintf('Constraining extrinsic signal to have mean 0!\n');           
            [pstar fstar eflag output] = fmincon(@(p)ll_func2(p,Basepars,Stimpars,Trainpars),pstar,[],[],A,0,[],[],[],gopts); 
        else
            [pstar fstar eflag output] = fminunc(@(p)ll_func2(p,Basepars,Stimpars,Trainpars),pstar,gopts);
        end
    end
      
    [fstar gstar] = ll_func2(pstar,Basepars,Stimpars,Trainpars);
    [Hstar] = [];%ll_hess_mult(Hinfo,eye(size(gstar,1)),Basepars,Stimpars,Trainpars);

% IF Raw filter and we have a norm constraint in place in Basepars    
elseif ( (strcmp(Basepars.hessmode ,'mult') && filtermode && (isfield(Basepars,'normconstr') && Basepars.normconstr)) )
    % Stimulus filter norm-constraints in nonseparable mode 
    gopts_constr = optimset(gopts,'Hessian','user-supplied','HessMult',@(x,lambda,Y) ll_hess_mult_normconst(x,lambda,Y,Basepars,Stimpars,Trainpars),'GradConstr','on','Algorithm','interior-point','SubproblemAlgorithm','cg');
    [pstar fstar] = fmincon(@(p)ll_func2(p,Basepars,Stimpars,Trainpars),p0,[],[],[],[],[],[],@(x) norm_constr(x,Basepars,Stimpars,Trainpars),gopts_constr);    
    [fstar gstar] = ll_func2(pstar,Basepars,Stimpars,Trainpars);
    Hstar = [];
% Some cateory called seperabel -raw..   not entirely sure    
elseif (strcmp(Basepars.k_filtermode,'sep_raw') && (isfield(Basepars,'normconstr') && Basepars.normconstr))
    % Spatial filter norm constraints in separable mode    
    fprintf('CONSTRAINING NORMS OF SPATIAL FILTERS....\n');
    gopts_constr = optimset(gopts,'Hessian','user-supplied','HessFcn',@(x,lambda) hessian_sep_normconstr(x,lambda,Basepars,Stimpars,Trainpars),'GradConstr','on','Algorithm','interior-point','SubproblemAlgorithm','cg');
    [pstar fstar] = fmincon(@(p)ll_func2(p,Basepars,Stimpars,Trainpars),p02,[],[],[],[],[],[],@(x) norm_constr(x,Basepars,Stimpars,Trainpars),gopts_constr);    
    [fstar gstar Hstar] = ll_func2(pstar,Basepars,Stimpars,Trainpars);
% Rank 2!!!  this is where we will end up

                %check_gradhess(0.1.*randn(size(p0)),@(p)
        %ll_func2(p,Basepars,Stimpars,Trainpars), []);
        % Put a norm 1 constraint on the spatial filters
        if 0
        Basepars.normconstr_p = 1;
        Basepars.normconstr = 1;
        nc = @(x) norm_constr(x,Basepars,Stimpars,Trainpars);
        gopts_constr = optimset(gopts,'Hessian','user-supplied','HessFcn',@(x,lambda) hessian_sep_normconstr(x,lambda,Basepars,Stimpars,Trainpars),'GradConstr','on','Algorithm','interior-point','SubproblemAlgorithm','cg');
        pstar = normalize_sep_filts(pstar,Basepars,size(Trainpars.D,2),Basepars.normconstr_p);
        fprintf('Constraining Lp norm of spatial filters to be less than 1!\n');
        end

        % Linear bounds on the spatial filers
        if 0
        %lb = -inf.*ones(size(pstar));
        %ub = inf.*ones(size(pstar));
        %spaceidx = get_pars_idx(Basepars,1,size(Trainpars.D,2),'kspace');
        %lb(spaceidx) = -1;
        %ub(spaceidx) = 1;
        %pstar = bound_sep_filts(pstar,Basepars,size(Trainpars.D,2));
        end

        if 0
            % Fix the shape of the filter and optimize
            K = get_model_filters(pstar,1,Basepars,1);
            pstarfixfilt = [pstar(get_pars_idx(Basepars,1,size(Trainpars.D,2),'ext'));
                pstar(get_pars_idx(Basepars,1,size(Trainpars.D,2),'b'))
                1;
                pstar(get_pars_idx(Basepars,1,size(Trainpars.D,2),'ps'))];
            Basepars2 = Basepars; Basepars2.filtermode = 'fixfilt';
            Basepars2.K{1} = K;
            for j=1:Basepars.Nneurons
                Trainpars.kx = Stimpars.dt .* sum(fastconv(Stimpars.movie_ROI(Basepars.crop_idx(:,j),:),Basepars2.K{j},length(Basepars.crop_idx(:,j)),size(Stimpars.movie_ROI,2)))';
            end
            pfcn2 = setup_plotfcns(Basepars2.filtermode,Basepars2,Stimpars,Trainpars);
            gopts2 = optimset(gopts,'PlotFcns',pfcn2);
            [pstarfixfilt] = fminunc(@(p) ll_func2_fixfilt(pstarfixfilt,Basepars2,Stimpars,Trainpars), pstarfixfilt, gopts2);
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%   HERE IS FROM FILTER STIMULUS TRAIN
%%%%%%%%%%%%%%%%%%%%%%%%%%   %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%   %%%%%%%%%%%%%%%%%%%%%%