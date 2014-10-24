function doPariwise(celldat,train,test)

fit.g.type = 'exp';
fit.g.p = 0;
if celldat.polarity == -1
  fit.f.type = 'rectLinearNeg';
  fit.f.p = 0;
else
  fit.f.type = 'rectLinear';
  fit.f.p = 0;
end
v = getPairwiseIntAll(train,test,fit);
saveVmat(celldat,fit,v,'vMats');