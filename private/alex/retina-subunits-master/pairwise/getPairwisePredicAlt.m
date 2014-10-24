function out = getPairwisePredicAlt(prs,sig1,sig2,fit)

constants = prs(1);
weights = prs(2:5);

term1 = evalNonLin(weights(1)*sig1,fit.f);
term2 = evalNonLin(weights(2)*sig2,fit.f);
term3 = evalNonLin(weights(3)*sig1 + weights(4)*sig2,fit.f);
term4 = (constants(1));
out = evalNonLin(term1 + term2 + term3 + term4,fit.g);

etol = 10^-6;
iiz = out <= etol;
out(iiz) = etol;