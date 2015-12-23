function cg_plotpars(x,basepars,stimpars,trainpars)

nrows = basepars.Nneurons;
ncols = 6;
neff = size(trainpars.D,2);
npars_perneuron = 1+basepars.n*basepars.Mk+basepars.nofilters_postspike+(neff-1)*basepars.nofilters_coupling;
[logprob cifs kx] = train_ll3(x,basepars,stimpars,trainpars);
[g] = train_ll_grad_hessmult(x,basepars,stimpars,trainpars,cifs);
for j=1:basepars.Nneurons

    idx = (j-1)*npars_perneuron+1:j*npars_perneuron;
    
    subplot(nrows,ncols,(j-1)*ncols+1);
    cla;
    imagesc(reshape(x(idx(2:1+basepars.n*basepars.Mk)),basepars.n,basepars.Mk)), axis tight, axis image;

    subplot(nrows,ncols,(j-1)*ncols+2);
    xaxis = (0:trainpars.dt:(basepars.Mhist-1)*trainpars.dt);
    plot(xaxis,basepars.postspike_basis*x(idx(1+basepars.n*basepars.Mk+1:1+basepars.n*basepars.Mk+basepars.nofilters_postspike)),'b.-');
    set(gca,'XScale','log');
    if (size(trainpars.D,2) > 1)
        subplot(nrows,ncols,(j-1)*ncols+3);
        plot(xaxis,basepars.coupling_basis*reshape(x(idx(1+basepars.n*basepars.Mk+basepars.nofilters_postspike+1:end)),basepars.nofilters_coupling,neff-1));
        set(gca,'XScale','log');
    end
    
    subplot(nrows,ncols,(j-1)*ncols+4);
    
    plot_gof_wrapper(x(idx),basepars,stimpars,trainpars);

    subplot(nrows,ncols,(j-1)*ncols+5);
    optimValues.gradient = g;
    plot_gradient(x(idx),optimValues,basepars,stimpars,trainpars);
    
    subplot(nrows,ncols,(j-1)*ncols+6);
    plot_hessian(x(idx),basepars,stimpars,trainpars);
    
    
    
end
pause(0.1);

