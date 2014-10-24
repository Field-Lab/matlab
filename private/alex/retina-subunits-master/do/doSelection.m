function doSelection(dat)

%%
fid = fopen(fullfile(dat.figurePath,'pairwise','selections.txt'),'a');
fprintf(fid,'%g\r',dat.rgcId);

%%
v = loadVmat(dat,'vMats');
[vals inds] = sort(v(:),'descend');
[cone1 cone2] = ind2sub(size(v),inds);

coneWeights = abs(dat.staFit);
coneWeights = coneWeights - min(coneWeights(:));
coneWeights = coneWeights/max(coneWeights(:));
thresh = prctile(coneWeights,60);
[~,strong_inds] = sort(coneWeights,'descend');
%%
for i=1:min(length(strong_inds),6)
    fprintf(fid,'strong cone is #%03g\r',dat.coneIds(strong_inds(i)));
end

for i=1:min(length(cone1),10)
    bestLin1 = dat.coneIds(cone1(i));
    bestLin2 = dat.coneIds(cone2(i));
    if coneWeights(cone1(i)) > thresh & coneWeights(cone2(i)) > thresh & vals(i) > 5
        fprintf(fid,'good linearity pair is #%03g and #%03g, index is %g\r',bestLin1,bestLin2,vals(i));
    end
end
