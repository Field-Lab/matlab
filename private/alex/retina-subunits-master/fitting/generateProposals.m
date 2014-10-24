function proposals = generateProposals(fit)

vals = nchoosek([1:fit.S],2);
nProposals = size(vals,1);

for in=1:nProposals
    
    indsMerge = [vals(in,1) vals(in,2)];
    old = ones(1,fit.S);
    old(indsMerge) = 0;
    indsOld = find(old);
    
    tmpFit.B_s = fit.B_s(indsOld);
    tmpFit.J_s = fit.J_s(indsOld);
    tmpFit.I_sc = fit.I_sc(indsOld,:);
    tmpFit.A_sc = fit.A_sc(indsOld,:);
    
    tmpFit.B_s = [tmpFit.B_s mean(fit.B_s(indsMerge))];
    tmpFit.J_s = [tmpFit.J_s 1];
    tmpFit.I_sc = [tmpFit.I_sc; sum(fit.I_sc(indsMerge,:))];
    tmpFit.A_sc = [tmpFit.A_sc; sum(fit.A_sc(indsMerge,:))];
    
    tmpFit.A_sc = normMat(tmpFit.A_sc,'L1',2);
    
    tmpFit.S = size(tmpFit.I_sc,1);
    tmpFit.C = size(tmpFit.I_sc,2);
    tmpFit.f = fit.f;
    tmpFit.g = fit.g;
    tmpFit.error = fit.error;
    tmpFit.constraints.B_s = fit.constraints.B_s;
    tmpFit.constraints.A_sc = fit.constraints.A_sc;
    
    proposals(in).fit = tmpFit;
    proposals(in).change.type = 'merge';
    proposals(in).change.val = indsMerge;
 
end