function idx = get_train_idx(basepars,res)

if (~isfield(basepars,'leaveout_idx'))
    idx = 1:basepars.maxt;
else
    idx = setdiff(1:basepars.maxt,basepars.leaveout_idx);
end
idx = idx(:);

if (nargin > 1 && strcmp(res,'spike'))
    idx = reshape(repmat(basepars.fac*(idx'-1),basepars.fac,1) + repmat((1:basepars.fac)',1,length(idx)),length(idx)*basepars.fac,1);
end