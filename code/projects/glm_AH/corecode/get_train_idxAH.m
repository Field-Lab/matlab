function idx = get_train_idxAH(Basepars,res)

if (~isfield(Basepars,'leaveout_idx'))
    idx = 1:Basepars.maxt;
else
    idx = setdiff(1:Basepars.maxt,Basepars.leaveout_idx);
end
idx = idx(:);

if (nargin > 1 && strcmp(res,'spike'))
    idx = reshape(repmat(Basepars.spikebins_perstimframe*(idx'-1),Basepars.spikebins_perstimframe,1)...
        + repmat((1:Basepars.spikebins_perstimframe)',1,length(idx)),length(idx)*Basepars.spikebins_perstimframe,1);
end