function out = getPairwisePredicNull(prs,sig1,sig2,fit)

constants = prs(1);
weights = prs(2:3);

term1 = evalNonLin(weights(1)*sig1,fit.f);
term2 = evalNonLin(weights(2)*sig2,fit.f);
term3 = (constants(1));

out = evalNonLin(term1 + term2 + term3,fit.g);

etol = 10^-6;
iiz = out <= etol;
out(iiz) = etol;