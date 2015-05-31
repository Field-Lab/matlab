function [fitGMLM_f,output]= fitGMLM_full(fitGMLM,ss,mov_use)

%% Make finer spikes

stim_length=size(mov_use,2)/120;
% idx=1:stim_length*120; 
% % spikes= idx(binnedspikes==1)/120; % this part of code is wrong! obviously
% %cannot go from lower resolution to higher resolution !
%         spksGen = zeros(stim_length*1200,1);
%         for ispike=1:length(spikes)
%             spksGen(floor(spikes(ispike)*1200)+1)=1;
%         end
%         spksGen2 = spksGen(1:stim_length*1200);
   
        
       spksGen=ss;    
%% Filter spikes and make history input 
% histBasLen=8;
% pars.ncols=histBasLen;
% pars.hpeaks=[2,30];
% pars.b=1;
% 
% dt=1;
% [iht, ihbas, ihbasis] = makeBasis_PostSpike(pars,dt);
% figure; 
% plot(iht/1200,ihbasis);
% title('Post Spike filter basis - Ask around for best thing');
% 
% 
% spks2=zeros(length(spksGen)+length(iht)-1,1);
% spks2(length(iht):length(spksGen)+length(iht)-1,1)=spksGen; % Append zeros before the movie
% 
% histInp=zeros(pars.ncols,length(spksGen));
% for ibasis=1:pars.ncols
% histInp (ibasis,:) = conv(spks2,ihbas(:,ibasis),'valid');
% end


% Perhaps we should combine this! With convolving with spikes !
    GLMPars = GLMParams;
    bin_size      = 1/1200;
    basis_params  = GLMPars.spikefilters.ps;
    ihbasis      = prep_spikefilterbasisGP(basis_params,bin_size);
    histBasLen=20;
    
    t_bin        = 1/1200;
    bins=stim_length*1200;
%  home_sptimes = spikes;
%  home_spbins  = ceil(home_sptimes / t_bin);
%  home_spbins = home_spbins(find(home_spbins < bins) );
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

    if(isfield(fitGMLM.Linear,'bias')) % if bias is fitted separately, then fit GMLM with added term.
        filteredStimDim = filteredStimDim+1;
        mov_filtered = [mov_filtered',ones(size(mov_filtered,2),1)]';

        for ifilter=1:nFrontEnds
            filters{ifilter} = [filters{ifilter};fitGMLM.Linear.bias{ifilter}];
        end

    end

% Approximate h .. based on regression .. 
h=0.01*ones( histBasLen,1);
% idxt=1:length(binnedResponses_global);
% tsp = idxt(binnedResponses_global>0);
% C= sum(hInp(:,tsp).*repmat(binnedResponses_global(tsp)',[histBasLen,1]),2)/(1/1200);
% 
% lam=exp(filters{1}'*mov_filtered);
% for ifilter=2:nFrontEnds
% lam = lam+ exp(filters{ifilter}'*mov_filtered);
% end
% lam = lam+mu;
% 
% B = sum(repmat(lam,[histBasLen,1]).*histInp,2);
% scale = ((B'*B)^(-1)) *B'*C;
% scale = log(scale/size(mov_filtered,2));
% h  = mean(hInp(:,tsp),2);
% h=h*scale/(norm(h)*sqrt(sum(hInp(:).^2)/size(mov_filtered,2)));
% %h=randn(histBasLen,1)/10;
figure;plot(ihbasis*h);title('Initial h');
pause(1);


    initialFilters=zeros(nFrontEnds*filteredStimDim +histBasLen+1,1);
    initialFilters(end)=mu;
    for ifilter=1:nFrontEnds
    initialFilters((ifilter-1)*filteredStimDim +1 : ifilter*filteredStimDim,1) = filters{ifilter};
    end
    initialFilters(ifilter*filteredStimDim+1:ifilter*filteredStimDim + histBasLen) = h;
    
%% Fit it.

optim_struct = optimset(...
   'derivativecheck','off',...
   'diagnostics','off',...  % 
   'display','iter',...  %'iter-detailed',... 
   'funvalcheck','off',... % don't turn this on for 'raw' condition (edoi).
   'GradObj','on',...
   'largescale','on',...
   'Hessian','off');
 [x,fval,exitflag,output,grad,hessian]  = fminunc(@GMLM_full_fcn,initialFilters,optim_struct);
 mu = x(end)

filters=cell(nFrontEnds,1);
for ifilter=1:nFrontEnds
filters{ifilter}=x((ifilter-1)*filteredStimDim+1:ifilter*filteredStimDim);
end


    histBas = x(nFrontEnds*filteredStimDim+1:nFrontEnds*filteredStimDim + histBasLen);
    histexpanded = ihbasis * histBas;

fitGMLM_f.Linear.filter=filters;
fitGMLM_f.mu=mu;
fitGMLM_f.hist.hBas= histBas;
fitGMLM_f.hist.hexpanded = histexpanded;
fitGMLM_f.hist.hbasis = ihbasis;


end