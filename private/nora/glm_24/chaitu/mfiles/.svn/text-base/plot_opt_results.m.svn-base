function [] = plot_opt_results(Pstar,p0, basepars,stimpars,trainpars,filtermode,f,g,Curvatures)

% Plot useful properties of optimization results

% 1st row: linear filter (spatial projection), linear filter (temporal properties), postspike filter
% 2nd row: sens

pstar = Pstar(:,1);

if (~isfield('psbasisGrad','trainpars'))
    trainpars.psbasisGrad = grad_basis(trainpars.D,basepars.postspike_basis);
end

if (strcmp(filtermode,'nonsep'))
    trainpars.lgrad = nonsep_lgrad(basepars,stimpars,trainpars);
end

ncols = 0;
npars = get_npars(basepars,size(trainpars.D,2));
switch(filtermode)
    case 'nonsep'
        nrows = 2*basepars.Nneurons;
        ncols = 5;
    case 'sep_raw'
        nrows = 3;
        ncols = 4;
        nspace = basepars.n;
        ntime = basepars.Mk;
    case 'sep_basis'
        nrows = 3;
        ncols = 4;
        nspace = basepars.nofilters_kspace;
        ntime = basepars.nofilters_ktime;
end

for j=1:basepars.Nneurons
    
    
    offset = (j-1)*2*ncols+1;
    
    
    [K PS C] = get_model_filters(pstar,j,basepars,size(trainpars.D,2));
    [K0 PS0 C0] = get_model_filters(p0,j,basepars,size(trainpars.D,2));
    xaxis = 0:stimpars.dt:(basepars.Mk-1)*stimpars.dt;    
    
    %K = reshape(basepars.ACmat*K(:),size(K));
    %K0 = reshape(basepars.ACmat*K0(:),size(K0));
    
    if (strcmp(filtermode,'nonsep'))
        
   
        1;
        
        if (basepars.n > 1)

            % Do the same for the optimized filter
            subplot(nrows,ncols,offset)
            imagesc(K), colorbar, xlabel(strcat('time (up to ',num2str((basepars.Mk-1)*stimpars.dt),' sec)')),ylabel('space'), title('Optimized linear filter in space-time');
            axis tight, axis image;    
            offset = offset + 1;            
            
            % Show Rank 1 approximation of the space-time stimulus filter
            [U S V] = svd(K);
            [spacefilt tempfilt] = get_spacetime_approx(K,1,0);
            
            
            %spacefilt = reshape(S(1,1).*U(:,1),basepars.box_m,basepars.box_n);
            %tempfilt = V(:,1);
            
            % There is indeterminacy in sign. Make sure it's positive when both are negative
            %if (sum(spacefilt(:))<0 && sum(tempfilt(:))<0)
            %    spacefilt = -spacefilt;
            %    tempfilt = -tempfilt;
            %end
            
            % Spatial filter
            subplot(nrows,ncols,offset)
            imagesc(spacefilt), colorbar, title(sprintf('Rank-1 spatial filter (%f percent)',S(1,1)./sum(diag(S))*100)), xlabel('window width'), ylabel('window height');
            axis image;
            offset = offset + 1;
            
            % Temporal filter
            subplot(nrows,ncols,offset)
            xaxis = 0:stimpars.dt:(basepars.Mk-1)*stimpars.dt;
            plot(xaxis,tempfilt), title('Rank-1 temporal filter'), xlabel('time preceding stimulus (s)'), ylabel('filter value');
            hold on, plot(xaxis,0,'--');
            offset = offset + 1;
            
            
        else
            
            subplot(nrows,ncols,offset), plot(0:stimpars.dt:(basepars.Mk-1)*stimpars.dt,K), title('Stimulus filter'), xlabel('time preceding stimulus (s)'), ylabel('filter value');
            offset = offset + 1;
            
        end
        
        
    else
        
        % Filters are separable - retrieve them
        [space1 temp1 space2 temp2] = get_sep_filters(p0,nspace,ntime,basepars.kspace_basis,basepars.ktime_basis);
        [space1star temp1star space2star temp2star] = get_sep_filters(pstar,nspace,ntime,basepars.kspace_basis,basepars.ktime_basis);
        % Plot the initial spatial filters
        subplot(nrows,ncols,offset);
        
        imagesc(reshape(space1,basepars.box_m,basepars.box_n)), title('Initial spatial filter 1'), xlabel('width'), ylabel('height'), colorbar;
        offset = offset + 1;
        subplot(nrows,ncols,offset);
        imagesc(reshape(space2,basepars.box_m,basepars.box_n)), title('Initial spatial filter 2'), xlabel('width'), ylabel('height'), colorbar;
        offset = offset + 1;
        % Plot the optimized spatial filters
        subplot(nrows,ncols,offset);
        imagesc(reshape(space1star,basepars.box_m,basepars.box_n)), title('Optimized spatial filter 1'), xlabel('width'), ylabel('height'), colorbar;
        offset = offset + 1;
        subplot(nrows,ncols,offset);
        imagesc(reshape(space2star,basepars.box_m,basepars.box_n)), title('Optimized spatial filter 2'), xlabel('width'), ylabel('height'), colorbar;
        offset = offset + 1;
        
        subplot(nrows,ncols,offset);
        plot(xaxis,[temp1 temp2],'.-'), title('Initial Temporal filters'), xlabel('time (s)'), ylabel('filter value');
        if (isfield(basepars,'ktime_basis') && ~isempty(basepars.ktime_basis))
            scl = max(abs([temp1;temp2]))/(2*max(abs(basepars.ktime_basis(:))));
            hold on, plot(xaxis,scl.*basepars.ktime_basis,'--');
        end
        offset = offset + 1;
        subplot(nrows,ncols,offset);
        plot(xaxis, [temp1star temp2star],'.-'), title('Optimized Temporal filters'), xlabel('time (s)'), ylabel('filter value');
        if (isfield(basepars,'ktime_basis') && ~isempty(basepars.ktime_basis))
            scl = max(abs([temp1star;temp2star]))/(2*max(abs(basepars.ktime_basis(:))));
            hold on, plot(xaxis,scl.*basepars.ktime_basis,'--');
        end        
        offset = offset + 1;
        
    end
        
    
    % Plot the history filter(s)
    psaxis = 0:trainpars.dt:(basepars.Mhist-1)*trainpars.dt;
    
    subplot(nrows,ncols,offset);
    plot(psaxis,PS,'.-'), hold on, plot(psaxis,PS0,'r--');%, legend('optimized','initial');
    set(gca,'XScale','log');
    title('Postspike filter'), xlabel('time after spike (s)'), ylabel('filter value');
    % Superimpose the basis fns
    scl = max(abs(PS))/max(abs(basepars.postspike_basis(:)));
    hold on, plot(psaxis,basepars.postspike_basis*scl,'--');
    
    offset = offset + 1;
    
    
    if (basepars.nofilters_coupling > 0 && size(trainpars.D,2) > 1)
        cpaxis = 0:trainpars.dt:(basepars.Mcoup-1)*trainpars.dt;
        subplot(nrows,ncols,offset)
        plot(cpaxis,C,'.-');
        set(gca,'XScale','log');
        title('Coupling filters'), xlabel('time after spike (s)'), ylabel('filter value');
        scl = max(abs(C(:)))/max(abs(basepars.postspike_basis(:)));
        hold on, plot(cpaxis,basepars.coupling_basis*scl,'--');
    
        offset = offset + 1;
    end
    
  
    % Plot some sensitivity values
    % trainpars.hst = hess_spacetimeconst(stimpars.x,basepars.ktime_basis,stimpars.dt);

    if (size(Pstar,2) > 1) % Plot the bootstrapping result
        subplot(nrows,ncols,offset);
        %size(Pstar)
        
        ps_mean = mean(Pstar(end-basepars.nofilters_postspike+1:end,2:end),2);
        ps_std = std(Pstar(end-basepars.nofilters_postspike+1:end,2:end),0,2);
        bar(ps_mean);
        hold on, errorbar(ps_mean,ps_std,'r.');
        title('Postspike basis coeff mean/std. dev. (bootstrapping)'), xlabel('ps filter number');
        offset = offset + 1;
    end    
    
    fprintf('Norm of the gradient is %f\n',norm(g));
    
    
    
    % How well-constrained are the parameters? Check the eigenvectors of the
    % Hessian to see curvature at the solution
    
    idx = (j-1)*npars+1:j*npars;
    
    subplot(nrows,ncols,offset);
    optimValues.gradient = g;
    plot_gradient(pstar(:,1),optimValues,basepars,stimpars,trainpars);
    offset = offset + 1;
    title('First-order perturbation analysis');
    switch(basepars.filtermode)
        case 'nonsep'
            xlabel(sprintf('base rate, mean stim filt, %d ps filts, %d cp filts',basepars.nofilters_postspike,basepars.nofilters_coupling));
        case 'sep_raw'
            xlabel(sprintf('base rate,sp1,t1,sp2,t2,%d ps filts, %d cp filts',basepars.nofilters_postspike,basepars.nofilters_coupling));
        case 'sep_basis'
            xlabel(sprintf('base rate,%d sp1, %d t1, %d sp2, %d t2,%d ps filts, %d cp filts',basepars.nofilters_kspace,basepars.nofilters_ktime,basepars.nofilters_kspace,basepars.nofilters_ktime,basepars.nofilters_postspike,basepars.nofilters_coupling));
    end
    ylabel('log|Df|');        

    subplot(nrows,ncols,offset);
    plot_hessian(pstar(:,1),basepars,stimpars,trainpars);
    title('Second-order perturbation analysis');
    switch(basepars.filtermode)
        case 'nonsep'
            xlabel(sprintf('base rate, mean stim filt, %d ps filts, %d cp filts',basepars.nofilters_postspike,basepars.nofilters_coupling));
        case 'sep_raw'
            xlabel(sprintf('base rate,sp1,t1,sp2,t2,%d ps filts, %d cp filts',basepars.nofilters_postspike,basepars.nofilters_coupling));
        case 'sep_basis'
            xlabel(sprintf('base rate,%d sp1, %d t1, %d sp2, %d t2,%d ps filts, %d cp filts',basepars.nofilters_kspace,basepars.nofilters_ktime,basepars.nofilters_kspace,basepars.nofilters_ktime,basepars.nofilters_postspike,basepars.nofilters_coupling));
    end   
    ylabel('log|Df|'); 
    offset = offset + 1;
    
    
    
    %[f g H] = ll_func2(Pstar(:,1),basepars,stimpars,trainpars);
    %[V D] = eig(H);
    %[sorted sorted_idx] = sort(abs(diag(D)),'ascend');
    %d = diag(D);
    %d = d(sorted_idx);
    
    % Plot the spectrum - this is not that informative
    %subplot(nrows,ncols,offset+2);
    %plot(log(abs(d))), title('log spectrum of Hstar'), xlabel('eigenvalue number'), ylabel('log(|eval|)');
    
    % Plot the eigenvectors
    %subplot(nrows,ncols,offset+4);
    %showIm(V(:,sorted_idx)), colormap(jet), colorbar, xlabel('eigenvector (incr. eval)'), ylabel('component val'), title('Eigendecomposition of Hstar');
    
    % Take the first r eigenvectors and plot their component values
    %r = 3;
    %subplot(nrows,ncols,offset+4);
    %plot(V(:,sorted_idx(1:r))), title(strcat('First ',num2str(r),' evecs')), xlabel('parameter index'), ylabel('evector val');
    
    % Perturbation analysis
    %del = 10^(-10); % perturbation magnitude
    
    % the ith component of Df is the approximate perturbation in the function if we
    % perturb by del in the ith direction, according to quadratic approx at the point
    %Df = Curvatures(:,1);
    %Df = del.*g + del^2 .* diag(H);
    %subplot(nrows,ncols,offset+5);
    %plot(log(abs(Df))), title('Perturbation analysis'), xlabel('parameter index'), ylabel('log(|DF|)');
    %hold on, plot(log(abs(g)),'r-');
    %legend('2nd order','1st order');
    
    % Plot goodness of fit
    subplot(nrows,ncols,offset);
    [lg_p cifs kx] = train_ll3(pstar,basepars,stimpars,trainpars);
    [uISIs uRISIs uLISIs KS] = goodness_of_fit(cifs(:,j),reprows(kx(:,j),basepars.fac),pstar((j-1)*npars+1),trainpars.dt,trainpars.D(:,trainpars.baseneuron_idx(j)),1);
    yaxis = linspace(0,1,length(uRISIs{1}))';
    plot(yaxis,yaxis,'b-');
    hold on, plot(yaxis,uISIs{1},'g--');
    %hold on, plot(yaxis,uLISIs{1},'r--');
    hold on, plot(yaxis,uRISIs{1},'r--');
    title(sprintf('CDF of rescaled ISIs. Nspikes=%d',length(uRISIs{1}))), xlabel(sprintf('K-S test: %d',KS(1)));
    ci = 1.36/sqrt(length(uRISIs{1})); % 95 % conf. interval
    tidx = ((yaxis + ci) < 1);
    bidx = ((yaxis - ci) > 0);
    hold on, plot(yaxis(tidx),yaxis(tidx)+ci,'m-');
    hold on, plot(yaxis(bidx),yaxis(bidx)-ci,'m-');    
end