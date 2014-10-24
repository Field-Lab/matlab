function selectCones(dat,v,file)

% select cone pairs worth running based on a data structure
% and a matrix of linearity indices (assumes upper triangular)

%%
fid = fopen(file,'w');
% sort the linearity indices
[vals inds] = sort(v(:),'descend');
[cone1 cone2] = ind2sub(size(v),inds);
coneWeights = abs(dat.staFit);
coneWeights = coneWeights - min(coneWeights(:));
coneWeights = coneWeights/max(coneWeights(:));

fprintf('\r')

for i=1:3
bestLin1 = dat.coneIds(cone1(i));
bestLin2 = dat.coneIds(cone2(i));
fprintf(fid,'good linearity pair is #%03g and #%03g, index is %g\r',bestLin1,bestLin2,vals(i));
end
%%

% sort by average cone strength
fprintf('\r');
flag = 1;
diffThresh = 0.2;
wThresh = 0.3;
linThresh = 2;
[wvals winds] = sort((coneWeights(cone1)+coneWeights(cone2))/2,'descend');
for i=1:length(winds)
    if cone1(winds(i)) < cone2(winds(i))
        if wvals(i) > wThresh
            if abs(coneWeights(cone1(winds(i))) - coneWeights(cone2(winds(i)))) < diffThresh;
                if abs(vals(winds(i))) < linThresh;
                    fprintf(fid,'good non linearity pair is #%03g and #%03g, index is %g\r',dat.coneIds(cone1(winds(i))),dat.coneIds(cone2(winds(i))),vals(winds(i)));
                end
            end
        end
    end
end

fclose(fid);


