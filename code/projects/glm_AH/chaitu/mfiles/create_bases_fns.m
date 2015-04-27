function [basepars trainpars]  = create_bases_fns(basepars0,stimpars,trainpars0)

basepars = basepars0;
trainpars = trainpars0;

% Number of timebins used to store postspike/coupling filter
basepars.Mhist = floor((basepars.fac*basepars.Mk)*basepars.ps_fratio); % on fine timescale ~ 8/factor m(s)^-1
basepars.Mcoup = floor((basepars.fac*basepars.Mk)*basepars.cp_fratio); % on fine timescale ~ 8/factor m(s)^-1
basepars.Nneurons = length(cells_totrain); % number of neurons fit in each training iteration
trainpars.dt = stimpars.dt / basepars.fac; % time resolution of spikes

% Set up the basis functions to use for training - uses the fac,fratio
% arguments to set up the timescales etc.

% K temporal basis parameters
alpha_ktime = 0.002;
beta_ktime = basepars.Mk*stimpars.dt; % ending pt of k filters

%bstretch = 0.002;
bstretch = basepars.bstretch;
%bstretch = 0.25; % stretch factor inside the log - this is the DEFAULT VALUE

% Linear filters
if (strcmp(filtermode,'sep_basis'))
    if (~isfield(basepars,'ktime_basis') || isempty(basepars.ktime_basis))
        fprintf('Creating temporal k basis.\n');
1;
        basepars.ktime_basis = create_histbasis(alpha_ktime,beta_ktime,bstretch,basepars.nofilters_ktime,basepars.Mk,stimpars.dt,'L2',basepars.k_spacing);
    else
        fprintf('Using supplied temporal k basis.\n');
    end
    if (~isfield(basepars,'kspace_basis'))
        fprintf('Creating spatial k basis.\n');
        % K spatial basis parameters
        min_std = basepars.min_std;
        max_std = basepars.max_std;
        centers = basepars.kspace_centers;
        basepars.kspace_basis = construct_gaussian_bump_basis(centers,min_std,max_std,(basepars.nofilters_kspace-1)/size(centers,1),sqrt(basepars.n),0);
        % add a background level (flat level) basis
        basepars.kspace_basis = [basepars.kspace_basis ones(basepars.n,1)];
        
        % Check that its the right size
        if (size(basepars.kspace_basis,2) ~= basepars.nofilters_kspace)
            error('Please check the spatial basis parameters. Make sure to leave one more for the flat basis.\n');
        end
        
        % Normalize
        basepars.kspace_basis = basepars.kspace_basis./repmat(sqrt(sum(basepars.kspace_basis.^2)),size(basepars.kspace_basis,1),1);
        
    else
        fprintf('Using supplied spatial k basis.\n');
    end
    

else
        basepars.ktime_basis = [];
        basepars.kspace_basis = [];
end

% PS/CP basis parameters (WARNING: magic numbers!)
alpha_ps = basepars.alpha_ps; %0.0; % starting point for PS
alpha_cp = basepars.alpha_cp; %-0.01; % starting pt for CP
beta_ps = (basepars.Mhist-1)*trainpars.dt; % ending point for ps filters
beta_cp = (basepars.Mcoup-1)*trainpars.dt; % '' cp filters


        
% History filters (postspike and coupling)
% Postspike filters

if (~isfield(basepars,'ps_spacing'))
    basepars.ps_spacing = pi/2;
end

basepars.postspike_basis = create_histbasis(alpha_ps,beta_ps,bstretch,basepars.nofilters_postspike,basepars.Mhist,trainpars.dt,'L2',basepars.ps_spacing);
%basepars.postspike_basis = create_exponential_basis(linspace(0.001,1/2*(basepars.Mhist*trainpars.dt),basepars.nofilters_postspike),basepars.Mhist,trainpars.dt,'L2');
% Coupling filters
if (find(cells_touse ~= cells_totrain))
    basepars.coupling_basis = create_histbasis(alpha_cp,beta_cp,bstretch,basepars.nofilters_coupling,basepars.Mcoup,trainpars.dt,'L2',basepars.ps_spacing);
else
    basepars.coupling_basis = [];
end

plot_basisfns(basepars.ktime_basis, basepars.kspace_basis,basepars.postspike_basis, basepars.coupling_basis, stimpars.dt, trainpars.dt,'linear');

1;

if 0 % Orthogonalize bases!
    
    
    fprintf('ORTHOGONALIZING THE BASES SETS!\n');
    
    if (isfield(basepars,'postspike_basis') && ~isempty(basepars.postspike_basis))
        basepars.postspike_basis = orthogonalize_basis(basepars.postspike_basis);
    end
    if (isfield(basepars,'coupling_basis') && ~isempty(basepars.coupling_basis))
        basepars.coupling_basis = orthogonalize_basis(basepars.coupling_basis);
    end
    if (isfield(basepars,'kspace_basis') && ~isempty(basepars.kspace_basis))
        basepars.kspace_basis = orthogonalize_basis(basepars.kspace_basis);
    end
    if (isfield(basepars,'ktime_basis') && ~isempty(basepars.ktime_basis))
        basepars.ktime_basis = orthogonalize_basis(basepars.ktime_basis);
    end
end
    