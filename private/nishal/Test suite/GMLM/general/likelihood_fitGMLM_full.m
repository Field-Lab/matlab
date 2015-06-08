function [nlikelihood,lam]= likelihood_fitGMLM_full(fitGMLM,ss,mov_use)



stim_length=size(mov_use,2)/120;

        
       spksGen=ss;    
%% Filter spikes and make history input 
%

% Perhaps we should combine this! With convolving with spikes !
    GLMPars = GLMParams;
    bin_size      = 1/1200;
    basis_params  = GLMPars.spikefilters.ps;
    ihbasis      = prep_spikefilterbasisGP(basis_params,bin_size);
    histBasLen=20;
    
    t_bin        = 1/1200;
    bins=stim_length*1200;

idx=1:length(spksGen);
home_spbins = idx(spksGen>0);

    basis         = ihbasis';
   histInp       = prep_convolvespikes_basis(home_spbins,basis,bins);

%% initialize with fitGMLM filters
global_vars_GMLM_afterSTC
% Expand mov_use 10 times
mov_use_big = zeros(size(mov_use,1),10*size(mov_use,2));

for istimdim=1:size(mov_use,1)
xx=mov_use(istimdim,:);
xx=repmat(xx,[10,1]);
xx=xx(:)';
mov_use_big(istimdim,:)=xx;
end

mov_filtered=mov_use_big;
binnedResponses_global=spksGen;
hInp = histInp;

filteredStimDim = length(fitGMLM.Linear.filter{1});
nFrontEnds = length(fitGMLM.Linear.filter);

mu=fitGMLM.mu;
filters=fitGMLM.Linear.filter;

    if(filteredStimDim > size(mov_filtered,1)); % if bias is fitted separately, then fit GMLM with added term.
     
        mov_filtered = [mov_filtered',ones(size(mov_filtered,2),1)]';

    end

% Approximate h .. based on regression .. 
h=fitGMLM.hist.hBas;

figure;plot(ihbasis*h);title('h');
pause(1);


    initialFilters=zeros(nFrontEnds*filteredStimDim +histBasLen+1,1);
    initialFilters(end)=mu;
    for ifilter=1:nFrontEnds
    initialFilters((ifilter-1)*filteredStimDim +1 : ifilter*filteredStimDim,1) = filters{ifilter};
    end
    initialFilters(ifilter*filteredStimDim+1:ifilter*filteredStimDim + histBasLen) = h;
    
%% Fit it.
 [nlikelihood,Grad,lam] = GMLM_full_fcn_return_stuff(initialFilters);


       
end