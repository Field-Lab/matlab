% THIS IS ALMOST AN UNNECESARRY AND USELES IF WE ARE ASSUMING RANK-2
% AH starting to edit on 10-6
% Training the parameters of the glm to data given a fixed set of stimulus filters

% Parameters to train:
% 1. stimulus filter (space x time) for each neuron
% 2. the coupling filters h (for each pair of neurons) - this is represented as cooridinates wrt some finite basis (of exponentials).
%    This is a 3D matrix (basis-coordinate-vector x neuron 1 x neuron 2)
% 3. the base firing rate of each neuron. This is a vector (neurons x 1)

% The function to be optimized is the log-likelihood of the data:

% sum_i[ sum_{t_{i,j}}[ log(L_{i,t})] - sum_t[L_{i,t}dt] ]    %% what is i
% i in this case... repeats of how long or what?


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
    
% Trainpars - structure with the following fields:
    % dt - temporal resolution of spike train time bins
    % D - sparse logical matrix where ith column is the spike train of the ith neuron
    
% gopts - optimization options to pass to fminunc



% Returns:
% pstar - optimal solution returned by fminunc
% fstar, gstar,Hstar - function, gradient, and hessian value at the optimal point
% eflag,output - exit flag and output structure returned by fminunc

function [pstar fstar gstar Hstar eflag output] = train_glm2_AH(p0,Basepars,Stimpars,Trainpars,gopts,NITER)
% fprintf('check point @ train_glm2 (edoi)')
% keyboard

if (~exist('NITER','var'))
    NITER = 1;
end

%%% AH   rank2 os 0 ,, fulll rank or raw is a 1
filtermode = ( strcmp(Basepars.k_filtermode,'nonsep') || strcmp(Basepars.k_filtermode,'raw'));   %% make sure it s not something

if strcmp(Basepars.k_filtermode,'rk2')   %%%  AH for the rank2 we seemed to get sent here 
    
    outflag = 0;
    counter = 0;
    
    while (~outflag)
        pstar = p0;
        

    
        % Put a norm constraint on the PS filters if we wanted it
        if (isfield(Basepars,'psnormconstr') && Basepars.psnormconstr > 0)
            fprintf('Constraining the norm of the PS filter coefficients.\n');
            nc = @(x) ps_norm_constr(x,Basepars,Stimpars,Trainpars);
            gopts_constr = optimset(gopts,'Hessian','user-supplied','HessFcn',@(x,lambda) hessian_sep_psnormconstr(x,lambda,Basepars,Stimpars,Trainpars),'GradConstr','on','Algorithm','interior-point','SubproblemAlgorithm','cg');
        end

        %%%%%%%% WHERE THE PROGRAM ACTUALLY BEGINS FOR RANK 2, NO NORM CONSTRAINT   %%%%%%%%%
        
        %fprintf('Adjusting extrinsic signal if necessary.\n');
        %pstar = adjust_baserate_ext(pstar,Basepars,Stimpars,Trainpars);
        if (counter > 0)
            pstar = 0.1.*randn(size(pstar)); 
            gopts = optimset(gopts,'MaxIter',120); % allow more iterations
        end
        
        for j=1:NITER  %% NITER is by default 1  !!  
            %%% again so far this section is not applicable.. not sure what
            if (isfield(Basepars,'ext_timepts') && ~isempty(Basepars.ext_timepts))
                
                mtx = sum(interpmtx(Basepars.interpfilt,Basepars.ext_timepts,Basepars.maxt));
                A = zeros(1,length(pstar));
                A(get_pars_idx(Basepars,1,size(Trainpars.D,2),'ext')) = mtx;
                fprintf('Constraining extrinsic signal to have mean 0!\n');
                if (~exist('nc','var'))%
                    [pstar fstar eflag output] = fmincon(@(p) ll_func2AH(p,Basepars,Stimpars,Trainpars), pstar, [], [], A, 0, [], [], [], gopts);
                    % [pstar fstar eflag output] = fmincon(@(p) ll_func2(p,Basepars,Stimpars,Trainpars), pstar, [], [], A, 0, lb, ub, [], gopts);
                else
                    [pstar fstar eflag output] = fmincon(@(p) ll_func2AH(p,Basepars,Stimpars,Trainpars), pstar, [], [], A, 0, [], [], nc, gopts_constr);
                end 
            %%%%%%% this is where we get sent on the while loop after the first initial .. unconstrained opt  %%%%%  
            else
                if (~exist('nc','var'))
                    %[pstar fstar eflag output] = fmincon(@(p) ll_func2(p,Basepars,Stimpars,Trainpars), pstar, [], [], [], [], lb, ub, [], gopts);
                    fprintf('unconstrained optimization..')
                    [pstar fstar eflag output] = fminuncAH(@(p) ll_func2AH(p,Basepars,Stimpars,Trainpars),pstar,gopts);
                else
                    fprintf('constrained optimization..')
                    [pstar fstar eflag output] = fmincon(@(p) ll_func2AH(p,Basepars,Stimpars,Trainpars), pstar, [], [], [], [], [], [], nc, gopts_constr);
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
        
        %%%   do this before you get kicked out of loop for gstar/Hstar 
        [fstar gstar Hstar] = ll_func2AH(pstar,Basepars,Stimpars,Trainpars);
    end
end
