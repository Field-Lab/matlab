% Set up the plotting functions

function pfcn = setup_plotfcns(filtermode,basepars,stimpars,trainpars)

% ----------------------------
pext = @(x,optimValues,state) plot_ext_prog(x,1,basepars,stimpars,trainpars);
pf0 = @(x,optimValues,state) plot_fixfilt_prog(x,1,basepars,stimpars,trainpars);
pf1 = @(x,optimValues,state) plot_filter_prog(x,1,basepars,stimpars,trainpars);
psqf1 = @(x,optimValues,state) plot_sqfilter_prog(x,1,basepars,stimpars,trainpars);
psf1 = @(x,optimValues,state) plot_sep_filter_prog(x,1,basepars,stimpars,trainpars);
pssqf1 = @(x,optimValues,state) plot_sep_sqfilter_prog(x,1,basepars,stimpars,trainpars);
ppf1 = @(x,optimValues,state) plot_ps_filter(x,1,basepars,stimpars,trainpars);
pcf1 = @(x,optimValues,state) plot_cp_filter(x,1,basepars,stimpars,trainpars);
pgf1 = @(x,optimValues,state) plot_gof_wrapper(x,basepars,stimpars,trainpars,optimValues);
pgrd1 = @(x,optimValues,state) plot_gradient(x,optimValues,basepars,stimpars,trainpars);
phss1 = @(x,optimValues,state) plot_hessian(x,basepars,stimpars,trainpars);
%------------------------------

fprintf('Updating plot fcns...\n');

n = 1;
if (isfield(basepars,'ext_timepts') && ~isempty(basepars.ext_timepts))
    pfcn{n} = pext;
    n = n+1;
end


% Plot the current filter
switch(filtermode)
    case {'nonsep','raw'}
        pfcn{n} = pf1;
    case {'sep_raw','rk2'}
        pfcn{n} = psf1;
    case 'sep_basis'
        pfcn{n} = psf1;
    case 'fixfilt'
        pfcn{n} = pf0;
end

n = n+1;

if (isfield(basepars,'XsqK') && basepars.XsqK)
    switch(filtermode)
        case 'nonsep'
            pfcn{n} = psqf1;
        case 'sep_raw'
            pfcn{n} = pssqf1;
        case 'sep_basis'
            pfcn{n} = pssqf1;
    end
    n = n+1;
end

% Plot the postspike filter
pfcn{n} = ppf1;
% Plot the coupling filters, if any
pfcn{n+1} = pcf1;
% Plot the GOF
pfcn{n+2} = pgf1;
% Plot the Grad vals
pfcn{n+3} = pgrd1;
% Plot the hessian diag vals
pfcn{n+4} = phss1;