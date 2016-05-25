% Find a good initial point

% Returns
% sta - uncropped sta(s)
% rsta - cropped and rotated sta(s)
% cpidx - crop indices
% p0 - initial search point
% ACmat - stimulus autocorrelation matrix

%function [stas rstas cpidces p0 box_m box_n ACmat] = find_init_pt(sp_times,basepars,stimpars,trainpars,prelimIter,prelimFrames,rotateFlag)
function [stas rstas cpidces p0 box_m box_n ACmat par_ini] = find_init_pt(sp_times,basepars,stimpars,trainpars,prelimIter,prelimFrames,rotateFlag)

N = basepars.Nneurons; % number of neurons whose parameters are being trained
Neff = size(trainpars.D,2); % "" whose spikes are being used

filtermode = basepars.filtermode;

1;

% Get the STA and cropping indices for the N cells
[stas rstas cpidces box_m box_n ACmat] = crop_sta(sp_times,stimpars.x,basepars,stimpars,rotateFlag);
% edoi: check

%fprintf('sta is computed')
%keyboard

if (isempty(cpidces))
    p0  = 0;
    return;
end

% WARNING : from here on we will assume that all simultaneously trained
% cells will have the same spatial RF size (i.e. basepars.crop_idx will be
%a matrix of EQUALLY sized vectors

if (size(cpidces{1},1) ~= basepars.n)
    fprintf('Spatial dimension has been truncated to %d\n',size(cpidces{1},1));
    
    if (strcmp(basepars.filtermode,'sep_basis'))
        fprintf('Cannot do truncated RFs with separable basis filtermode. Skipping this cell.\n');
        cpidces = [];
        p0 = 0;
        return;
    end
    basepars.n = size(cpidces{1},1);
end

if (N > 1)
    fprintf('WARNING : assuming equally-sized RFs!\n');
end

cpidces = cell2mat(cpidces');
rstas = cell2mat(rstas');
stas = cell2mat(stas');


% Preliminary optimization to find initial point
basepars.crop_idx = cpidces;
basepars.maxt = prelimFrames;
stimpars.x = stimpars.x(:,1:prelimFrames,:);

1;

if 0 (strcmp(filtermode,'nonsep') && basepars.prelimIter > 0)
    
    % Assume raw separable filters in the preliminary mode: NOTE: changed (see next)
    basepars.filtermode = 'sep_raw';

    npars = 1+2*(basepars.n+basepars.Mk)+basepars.nofilters_postspike+(Neff-1)*basepars.nofilters_coupling;
    % Set the initial stim filter to the Rank-2 approx of the STA
    p0 = set_sta_init(rstas,basepars,npars,N,Neff);
    
    pfcn = setup_plotfcns(basepars.filtermode,basepars,stimpars,trainpars);
    gopts = optimset('GradObj','on','Hessian', 'on', 'TolFun',10^(-6),'MaxIter',prelimIter,'display','iter','PlotFcns',pfcn);
    switch basepars.prelim_mode
        case 'sep' % Separable mode 
            
            % Options
            % Run preliminary optimization
            fprintf('Running preliminary optimization (separable filters) for %d iterations.\n',prelimIter);
            
            1;
            
            [presult] = train_glm2(p0,basepars,stimpars,trainpars,gopts);
            nfullpars = (1+basepars.n*basepars.Mk+basepars.nofilters_postspike+(Neff-1)*basepars.nofilters_coupling);
            p0 = zeros(N*nfullpars,1);
            
            for j=1:N
                idx = (j-1)*npars+1:j*npars;
                fullidx = (j-1)*nfullpars+1:j*nfullpars;
                p0(fullidx(1)) = presult(idx(1));
            
                [s1 t1 s2 t2] = get_sep_filters(presult(idx),basepars.n,basepars.Mk);
            
                p0(fullidx(2:1+basepars.n*basepars.Mk)) = reshape((s1*t1') + (s2*t2'),basepars.n*basepars.Mk,1);
            
                p0(fullidx(1+basepars.n*basepars.Mk+1:end)) = presult(idx(1+2*(basepars.n+basepars.Mk)+1:end));
            end
            
        case 'mag' % Optimize the magnitude of the stimulus filter along with the history filters (currently unused)
             
            rsta_norm = rsta./norm(rsta); % set stim filter to normalized sta;
            % Optimize wrt stim filter magnitude, base rate, and postspike/coupling filters
            p0 = zeros(2+basepars.nofilters_postspike+(basepars.Nneurons-1)*basepars.nofilters_coupling,1);
            p0(2) = 1;
            % Initialize the trainpars structure
            if (~isempty(basepars.ktime_basis) && ~isfield(trainpars,'hst'))
                trainpars.hst = hess_spacetimeconst(stimpars.x,basepars.ktime_basis,stimpars.dt);
            end
            
            if (~isfield(trainpars,'psbasisGrad'))
                fprintf('Computing ps basis grad...');
                trainpars.psbasisGrad = grad_basis(trainpars.D,basepars.postspike_basis);
                fprintf('done.\n');
                %trainpars.psbasisGrad = trainpars.psbasisGrad(:,ps_
            end
            
            if (basepars.Nneurons > 1 && ~isfield(trainpars,'cpbasisGrad'))
                fprintf('Computing cp basis grad...');
                trainpars.cpbasisGrad = grad_basis(trainpars.D,basepars.coupling_basis);
                fprintf('done.\n');
            else
                trainpars.cpbasisGrad = [];
            end
            
            if (strcmp(basepars.filtermode,'nonsep') && ~isfield(trainpars,'lgrad'))
                fprintf('Computing lgrad for nonsep mode...');
                trainpars.lgrad = nonsep_lgrad(basepars,stimpars,trainpars);
                fprintf('done.\n');
            end
            
            % Filter the stimulus with the sta
            p0_full = magtofull(p0,basepars,rsta_norm(:));
            kx_sta = filterstimulus_train2_nonsep(p0_full,basepars,stimpars,trainpars);
            
            % Optimization use fminunc
            pfcn{1} = @(x,optimValues,state) plot_ps_filter(x,basepars,stimpars,trainpars);
            pfcn{2} = @(x,optimValues,state) plot_fixfilt_filter(x,rsta_norm,basepars,stimpars,trainpars);
            gopts = optimset(gopts,'HessMult',@(Hinfo,Y) ll_fixfilt_hess_mult(Hinfo,Y,basepars,stimpars,trainpars),'PlotFcns',pfcn);
            func = @(p)ll_func2_fixfilt_hessmult(p,rsta_norm,kx_sta,basepars,stimpars,trainpars);
            [p0] = fminunc(func,p0,gopts);
            fprintf('Completed optimization of postspike and base parameters\n');
            p0 = magtofull(p0,basepars,rsta_norm(:));
    end
    
    
end

if 1 %((strcmp(filtermode,'sep_raw') || strcmp(filtermode,'sep_basis')) && strcmp(basepars.prelim_mode,'nonsep')) % CHANGED: ED.
    bnsp = basepars;
    bnsp.filtermode = 'nonsep';
    npars_nonsep = get_npars(bnsp,Neff);    
    p0 = zeros(N*npars_nonsep,1);
    1;
    for j=1:N
        p0(get_pars_idx(bnsp,j,Neff,'k')) = rstas(:,j);
    end    
    
    % Turn off FREEZING CONSTRAINTS
    if (isfield(bnsp,'frozen_idx'))
        fprintf('DISABLING FROZEN PARAMETERS for preliminary optimization...\n');
        bnsp.frozen_idx = [];
        bnsp.frozen_vals = [];
    end
    
    
    % If any parameters are frozen at designated values, reset them to their
    % designated values
    1;
    frozenflag = isfield(bnsp,'frozen_idx') && isfield(bnsp,'frozen_vals');
    if (frozenflag)
        fprintf('Freezing %d parameters...\n',length(bnsp.frozen_vals));
        p0(bnsp.frozen_idx) = bnsp.frozen_vals;
    end

    % Run preliminary nonsep mode for a few iterations
    bnsp.maxNonsepIter = bnsp.prelimIter;
    
    1;
    pfcn = setup_plotfcns(bnsp.filtermode,bnsp,stimpars,trainpars);
    gopts = optimset('GradObj','on','Hessian', 'on', 'TolFun',10^(-6),'MaxIter',prelimIter,'display','iter','PlotFcns',pfcn);
    p0actual = 1e-5.*randn(size(p0));
    if(frozenflag)
        p0actual(bnsp.frozen_idx) = p0(bnsp.frozen_idx);
    end
    
    % TURN OFF REGULARIZATION IN PRELIM NONSEP MODE
    if (isfield(bnsp,'lambda_reg'))
        fprintf('DISABLING REGULARIZATION for preliminary optimzation...\n');
        bnsp.lambda_reg = 0;
    end
    
    % TURN OFF NORM CONSTRAINTS
    if (isfield(bnsp,'normconstr'))
        fprintf('DISABLING NORM CONSTRAINTS for preliminary optimization...\n');
        bnsp = rmfield(bnsp,'normconstr');
    end
    

    if (size(stimpars.x,3)>1)
        1;
        pnonsep = train_glm_nonsep_multi(p0actual,bnsp,stimpars,trainpars,gopts);
    else
        [pnonsep,~,~,par_ini] = train_glm_nonsep(p0actual,bnsp,stimpars,trainpars,gopts); % edoi, chk, 2012-01-08
    end
    
    %keyboard;
    
    npars_sep = get_npars(basepars,Neff);
    p0 = zeros(npars_sep,1);    
    
    %keyboard;
    
    for j=1:N
        
        
        if (isfield(basepars,'ext_timepts') && ~isempty(basepars.ext_timepts))
            p0(get_pars_idx(basepars,j,Neff,'ext')) = pnonsep(get_pars_idx(bnsp,j,Neff,'ext')); 
            fprintf('ind_init_pt check (edoi)')
            keyboard 
        end
        
        p0(get_pars_idx(basepars,j,Neff,'b')) = pnonsep(get_pars_idx(bnsp,j,Neff,'b'));
        
        
        switch(filtermode)
           case 'nonsep'
              p0(get_pars_idx(basepars,j,Neff,'k')) = get_model_filters(pnonsep,j,bnsp,Neff);
              
           case 'sep_raw'
              [space time] = get_spacetime_approx(get_model_filters(pnonsep,j,bnsp,Neff),2,0);
              p0(get_pars_idx(basepars,j,Neff,'k')) = [space(:,1);time(:,1);space(:,2);time(:,2)];  % copy spacetime filters
              
              if (isfield(basepars,'XsqK') && basepars.XsqK)
                 [blah blah blah Ksq] = get_model_filters(pnonsep,j,bnsp,Neff);
                 [space time] = get_spacetime_approx(Ksq,2,0);
                 p0(get_pars_idx(basepars,j,Neff,'ksq')) = [space(:,1);time(:,1);space(:,2);time(:,2)];  % copy spacetime filters
              end
              
              
              
           case 'sep_basis'
              [space time] = get_spacetime_approx(get_model_filters(pnonsep,j,bnsp,Neff),2,0);
              lambda = 0.0;
              
              p0(get_pars_idx(basepars,j,Neff,'k')) = [project_to_subspace(basepars.kspace_basis,space(:,1),lambda);...
                 project_to_subspace(basepars.ktime_basis,time(:,1),lambda);...
                 project_to_subspace(basepars.kspace_basis,space(:,2),lambda);...
                 project_to_subspace(basepars.ktime_basis,time(:,2),lambda)];
              
              if (isfield(basepars,'XsqK') && basepars.XsqK)
                 [blah blah blah Ksq] = get_model_filters(pnonsep,j,bnsp,Neff);
                 [space time] = get_spacetime_approx(Ksq,2,0);
                 p0(get_pars_idx(basepars,j,Neff,'ksq')) = [project_to_subspace(basepars.kspace_basis,space(:,1),lambda);...
                    project_to_subspace(basepars,ktime_basis,time(:,1),lambda);...
                    project_to_subspace(basepars.kspace_basis,space(:,2),lambda);...
                    project_to_subspace(basepars,ktime_basis,time(:,2),lambda)];
              end
              
                
        end
        p0(get_pars_idx(basepars,j,Neff,'ps')) = pnonsep(get_pars_idx(bnsp,j,Neff,'ps')); % copy ps filters
        if (Neff>1)
            p0(get_pars_idx(basepars,j,Neff,'cp')) = pnonsep(get_pars_idx(bnsp,j,Neff,'cp')); % copy cp filters
        end
    end
    
    %keyboard;
    
end
if (strcmp(filtermode,'nonsep') && basepars.prelimIter <= 0)
   fprintf('@find_init_pt.m, nonsep mode, no prelim iter.\n')
   %keyboard
   nfullpars = get_npars(basepars,Neff);
   p0 = zeros(N*nfullpars,1);
   for j=1:N
      p0(get_pars_idx(basepars,j,Neff,'k')) = rstas(:,j);
   end
   p0 = 1e-5.*randn(size(p0)).*0.01;
end
%keyboard