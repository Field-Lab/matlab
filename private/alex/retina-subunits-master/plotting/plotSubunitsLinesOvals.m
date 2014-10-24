function plotSubunitsLinesOvals(fit,dat,coneSz,subSz,coneEdge)

% subfunction for plotModel2d.m that plots a subunit model fit
% using lines and ovals
%
% probably want to modify this so that it only does the plotting,
% not the creation of the rectangle or the setting of the axes

colorSub = 1;

if nargin < 3
    coneSz = 1.4;
    subSz = 2;
    coneEdge = 1.25;
end

centroid = getCentroid(dat.locs_c);

% plot the lines
locs_s = normMat(fit.I_sc,'L1',2)*dat.locs_c;
for is=1:fit.S
    startPt = [locs_s(is,1) locs_s(is,2)];
    endPt = centroid;
    h = line([startPt(1) endPt(1)],[startPt(2) endPt(2)]);
    weight = fit.B_s(is)/norm(fit.B_s);
    set(h,'Color','w');
    set(h,'LineWidth',abs(weight)*coneEdge*8);
end

% plot the subunit ovals
hold on
for is=1:fit.S
    if sum(fit.I_sc(is,:)) > 1
       %mu = mean(dat.locs_c(logical(fit.I_sc(is,:)),:));
       %c = cov(dat.locs_c(logical(fit.I_sc(is,:)),:));
       %c = c + subSz*eye(2);
       %h = drawEllipseFilled(mu,c,[0.75 0.75 0.75]);
       h = drawConvexShape(dat.locs_c(logical(fit.I_sc(is,:)),:),[0.75 0.75 0.75],subSz*0.5);
       set(h,'EdgeColor','none');
    else
        circle(locs_s(is,1),locs_s(is,2),subSz,[0.75 0.75 0.75]','none');
    end
end

if colorSub
    baseClrs = cmap_hsv_dark(10);
    subClrs = zeros(fit.S,3);
    nonSingCounts = sum(fit.I_sc,2);
    for j=2:max(nonSingCounts)
      tmpInds = sum(fit.I_sc,2) == j;
      if ~isempty(tmpInds)
        if j > 10
            val = 10;
        else
            val = j;
        end
        subClrs(tmpInds,:) = repmat(baseClrs(val-1,:),sum(tmpInds),1);
      end
    end
end

% plot the cones
hold on
xPos = dat.locs_c(:,1);
yPos = dat.locs_c(:,2);
for is=1:fit.S
    if sum(fit.I_sc(is,:)) > 1
        if colorSub
            inds = find(fit.I_sc(is,:));
            weights = fit.A_sc(is,inds);
            clr = repmat(subClrs(is,:),[length(weights) 1]);
            clr_hsv = rgb2hsv(clr);
            clr_hsv(:,2) = clr_hsv(:,2).*weights(:)/max(weights(:));
            clr_hsv(:,3) = clr_hsv(:,3)+0.2*(1-weights(:)/max(weights(:)));
            clr = hsv2rgb(clr_hsv);
            for j=1:length(inds)
                circle(xPos(inds(j)),yPos(inds(j)),coneSz,...
                    [clr(j,:)],'k',coneEdge);
            end
        else
            inds = find(fit.I_sc(is,:));
            weights = fit.A_sc(is,inds);
            weights = weights/max(weights);
            weights = 0.75 + weights*0.25;
            for j=1:length(inds)
                circle(xPos(inds(j)),yPos(inds(j)),coneSz,...
                    [weights(j) weights(j) weights(j)],'k',coneEdge);
            end
        end
    else
        inds = find(fit.I_sc(is,:));
        circle(xPos(inds(1)),yPos(inds(1)),coneSz,'w','k',coneEdge);
    end
end