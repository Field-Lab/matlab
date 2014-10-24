function init_orig = fitGreedy(celldat,train,test,plotting)

% start with singletons, evaluate fit quality
% try all possible pairs
% if we don't improve, stop
% lock that pair, then try all possible new pairs / triplets
% need a function where i give it an architecture and it does the fitting

%%
if exist('plotting','var') && plotting == 1
  figure;
  set(gcf,'Position',[598 460 842 346]);
end
%% start with singletons
%fprintf('(fitGreedy) starting singleton fit\n');
init_orig = mkSubAssgn(celldat.nConesFit,'singletons');
[bestfit negloglik_orig r2] = fitABFG(celldat,train,test,init_orig);
negloglik_iter = negloglik_orig;
%% try pairs
flag = 1;
iIter = 2;
while flag
  combos = combntns(1:size(init_orig,1),2);
  negloglik_try = zeros(size(combos,1),1);
  r2_try = zeros(size(combos,1),1);
  
  %fprintf('(fitGreedy) trying mergers\n');
  parfor_progress(size(combos,1));
  parfor i=1:size(combos,1)
    %i
    parfor_progress;
    init_try{i} = mergeRows(init_orig,[combos(i,1) combos(i,2)]);
    [bestfit_try{i} negloglik_try(i) r2_try(i)] = fitABFG(celldat,train,test,init_try{i},bestfit);
  end
  parfor_progress(0);
  
  if sum(negloglik_try < negloglik_orig)
    [~,mnInd] = min(negloglik_try);
    negloglik_orig = negloglik_try(mnInd);
    init_orig = init_try{mnInd};
    negloglik_iter(iIter) = negloglik_orig;
    bestfit = bestfit_try{mnInd};
    if exist('plotting','var') && plotting == 1
      subplot(1,2,1);
      axis equal
      xlim(celldat.xRange);
      ylim(celldat.yRange);
      set(gca,'YDir','reverse');
      axis off
      rect1 = rectangle('Position',[celldat.xRange(1) celldat.yRange(1)...
        sum(celldat.xRange) sum(celldat.yRange)]);
      set(rect1,'FaceColor',[0.5 0.5 0.5]);
      set(rect1,'EdgeColor','none');
      plotSubunitsLinesOvals(bestfit,celldat);
      subplot(1,2,2);
      plot([1:iIter],negloglik_iter(1:end),'.-');
      drawnow;
    end
    fprintf('(fitGreedy) merging subunits %g and %g\n',combos(mnInd,1),combos(mnInd,2));
    iIter = iIter + 1;
  else % no improvement in likelihood
    flag = 0;
    %fprintf('(fitGreedy) greedy search finished\n');
  end
  if size(init_orig,1) == 1 % all cones in a subunit
    flag = 0;
    %fprintf('(fitGreedy) greedy search finished\n');
  end
  %%
  
end

