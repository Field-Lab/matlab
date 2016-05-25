function out = plot_gof_wrapper(x,basepars,stimpars,trainpars,optimValues)

out = 0;
train_idx = get_train_idx(basepars);
if (strcmp(basepars.filtermode,'fixfilt'))
   [logprob lcifs cifs kx] = train_ll3_fixfilt(x,basepars,stimpars,trainpars);
   % Invert the leaveout
   if (isfield(basepars,'leaveout_idx'))
      basepars.leaveout_idx = setdiff(1:basepars.maxt,basepars.leaveout_idx);
      [cv_logprob] = train_ll3_fixfilt(x,basepars,stimpars,trainpars);
   end
else
   [logprob lcifs cifs kx] = train_ll3(x,basepars,stimpars,trainpars);
   
   % Invert the leaveout
   if (isfield(basepars,'leaveout_idx'))
      basepars.leaveout_idx = setdiff(1:basepars.maxt,basepars.leaveout_idx);
      [cv_logprob] = train_ll3(x,basepars,stimpars,trainpars);
   end
   
end
[uISIs uRISIs uLISIs KS nsp] = goodness_of_fit(cifs,reprows(kx,basepars.fac),x(1),trainpars.dt,trainpars.D(:,trainpars.baseneuron_idx),1);

if (isfield(basepars,'leaveout_idx') && ~isempty(basepars.leaveout_idx))
   title(sprintf('lgp=%0.3f, cvlgp=%0.3f, ks=%0.3f/%0.3f',logprob/length(train_idx),cv_logprob/(basepars.maxt-length(train_idx)),KS(1),1.36./nsp(1)));
else
   title(sprintf('lgp=%0.3f, ks=%0.3f/%0.3f',logprob/length(train_idx),KS(1),1.36./nsp(1)));
end


global CVLGP;
global CVLGP_OLD;


%if (cv_logprob < CVLGP && CVLGP < CVLGP_OLD && optimValues.funccount > 5) % over fitting , quit!
if ~isempty(CVLGP) % added, 2011/12/22
   if (cv_logprob < CVLGP && optimValues.funccount > 5)%% over fitting , quit!
      out = 1;
   end
end

CVLGP_OLD = CVLGP;
CVLGP = cv_logprob;
