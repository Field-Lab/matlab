function fit = addNeighbors(dat,fit,mult)

% fit = addNeighbors(dat,fit,thresh)
%
% given a set of assignments from cones to subunits,
% and the locations of all the cones,
% augment each subunit so that it contains
% a set of additional cones within about one cone distance
% from each of the current cones

locs_s = normMat(fit.I_sc,'L1',2)*dat.locs_c;
distMat_sc = getDistances(locs_s,dat.locs_c);
distMat_cc = getDistances(dat.locs_c,dat.locs_c);
min_dist = min(triuVec(distMat_cc,1));

for is=1:size(fit.I_sc,1)
  %if sum(fit.I_sc(is,:))>1
  %maxRad = max(distMat_sc(is,:));
  %inds = distMat_sc(is,:)<thresh*maxRad;
  %fit.I_sc(is,inds) = 1;
  %end
  add = zeros(1,size(fit.I_sc,2));
  inds = find(fit.I_sc(is,:)>0);
  if length(inds)>1
    for ic=1:length(inds)
      add = add + double(distMat_cc(inds(ic),:) > 0 & distMat_cc(inds(ic),:) < min_dist*mult);
    end
    fit.I_sc(is,:) = or(fit.I_sc(is,:),add>0);
  end
end
