function init_orig = fitGreedyLocal(dat,train,test,plotting)

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
%fprintf('(fitGreedyLocal) starting singleton fit\n');
init_orig = mkSubAssgn(dat.nConesFit,'singletons');
[bestfit negloglik_orig r2] = fitABFG(dat,train,test,init_orig);
negloglik_iter = negloglik_orig;
%% try pairs
flag = 1;
iIter = 2;
while flag
  combos = [];
  N_nearest = 6;
  locs_s = normMat(init_orig,'L1',2)*dat.locs_c;
  for is = 1:size(init_orig,1)
    current_inds = is;
    new_inds = [1:is-1 is+1:size(init_orig,1)];
    current_center = locs_s(is,:);
    new_center = locs_s(new_inds,:);
    d = getDistances(current_center,new_center);
    [~,inds] = sort(d,'ascend');
    combos_tmp = zeros(min(length(inds),N_nearest),2);
    combos_tmp(:,1) = is;
    combos_tmp(:,2) = new_inds(inds(1:min(length(inds),N_nearest)));
    combos = [combos; combos_tmp];
  end
  combos = sort(combos,2);
  combos = unique(combos,'rows');
  
  %%
  %combos = combntns(1:size(init_orig,1),2);
  negloglik_try = zeros(size(combos,1),1);
  r2_try = zeros(size(combos,1),1);
  
  %fprintf('(fitGreedyLocal) trying mergers\n');
  parfor i=1:size(combos,1)
    init_try{i} = mergeRows(init_orig,[combos(i,1) combos(i,2)]);
    [bestfit_try{i} negloglik_try(i) r2_try(i)] = fitABFG(dat,train,test,init_try{i},bestfit);
  end
  %%
  if sum(negloglik_try < negloglik_orig)
    [~,mnInd] = min(negloglik_try);
    negloglik_orig = negloglik_try(mnInd);
    init_orig = init_try{mnInd};
    negloglik_iter(iIter) = negloglik_orig;
    bestfit = bestfit_try{mnInd};
    %fprintf('(fitGreedyLocal) merging subunits %g and %g\n',combos(mnInd,1),combos(mnInd,2));
    if exist('plotting','var') && plotting == 1
      subplot(1,2,1);
      axis equal
      xlim(dat.xRange);
      ylim(dat.yRange);
      set(gca,'YDir','reverse');
      axis off
      rect1 = rectangle('Position',[dat.xRange(1) dat.yRange(1) sum(dat.xRange) sum(dat.yRange)]);
      set(rect1,'FaceColor',[0.5 0.5 0.5]);
      set(rect1,'EdgeColor','none');
      plotSubunitsLinesOvals(bestfit,dat);
      subplot(1,2,2);
      plot([1:iIter],negloglik_iter(1:end),'.-');
      drawnow;
    end
    iIter = iIter + 1;
  else % no improvement in likelihood
    flag = 0;
    %fprintf('(fitGreedyLocal) greedy search finished\n');
  end
  if size(init_orig,1) == 1 % all cones in a subunit
    flag = 0;
    %fprintf('(fitGreedyLocal) greedy search finished\n');
  end
  
  %%
  
end

