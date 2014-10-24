function proposals = evalProposals(train,proposals)

nIter = 2;
nProposals = length(proposals);
lik = zeros(1,nProposals);

for in=1:nProposals
    % optimize the weights for the candidate structure
    for iIter = 1:nIter
        proposals(in).fit = fitA(train,proposals(in).fit,0);
        proposals(in).fit = fitB(train,proposals(in).fit,0);
    end
    lik(in) = evalFit(train,proposals(in).fit);
    proposals(in).lik = lik(in);
    plotModel(proposals(in).fit); drawnow;
end

