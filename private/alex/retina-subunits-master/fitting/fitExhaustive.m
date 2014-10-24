function init_orig = fitExhaustive(celldat,train,test,plotting)

%% start with singletons
init_orig = mkSubAssgn(celldat.nConesFit,'singletons');
[bestfit negloglik_orig r2] = fitABFG(celldat,train,test,init_orig);
negloglik_iter = negloglik_orig;
%% try all possible combinations up to order n
flag = 1;
order = 3;
combos = foo; % get all possible combinations up to 'order'

for i=1:size(combos,1)
  init_try{i} = mergeRows(init_orig,[combos(i,1) combos(i,2)]);
  [bestfit_try{i} negloglik_try(i) r2_try(i)] = fitABFG(celldat,train,test,init_try{i},bestfit);
end

if sum(negloglik_try < negloglik_orig)
    [~,mnInd] = min(negloglik_try);
    negloglik_orig = negloglik_try(mnInd);
    init_orig = init_try{mnInd};
    negloglik_iter(iIter) = negloglik_orig;
    bestfit = bestfit_try{mnInd};
end