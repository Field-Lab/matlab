function plotPairwiseLines(dat,v,withNum,colored)

coneSz = 1.4;
coneEdge = 1.25;

% plot cones, colored by sta if we have it
if colored
    scatter(dat.locs_c(:,1),dat.locs_c(:,2),200,dat.staFit,'filled');
    set(gca,'CLim',[-max(abs(dat.staFit)) max(abs(dat.staFit))])
else
    xPos = dat.locs_c(:,1);
    yPos = dat.locs_c(:,2);
    for i=1:size(dat.locs_c,1)
        circle(xPos(i),yPos(i),coneSz,'w','k',coneEdge);
    end
end

if isfield(dat,'coneTypesFit')
    for i=1:size(dat.locs_c,1)
        if ~strcmp(dat.coneTypesFit(i),'U');
            h = text(dat.locs_c(i,1),dat.locs_c(i,2),dat.coneTypesFit(i));
            set(h,'FontWeight','bold');
            set(h,'HorizontalAlignment','center');
            set(h,'VerticalAlignment','middle');
        end
    end
end

% add a line between each cone with width
% given by the degree of linearity
vScaled = v/max(abs(v(:)));
vScaledVec = triuVec(vScaled,1);
vVec = triuVec(v,1);
colInds = meshgrid([1:size(v,1)]);
rowInds = colInds';
colInds = triuVec(colInds,1);
rowInds = triuVec(rowInds,1);
[~,inds] = sort(vScaledVec,'ascend');
for i=1:length(vScaledVec);
    if vScaledVec(inds(i)) > 0.3
        x = [dat.locs_c(rowInds(inds(i)),1) dat.locs_c(colInds(inds(i)),1)];
        y = [dat.locs_c(rowInds(inds(i)),2) dat.locs_c(colInds(inds(i)),2)];
        h = line(x,y);
        set(h,'Color','k');
        set(h,'LineWidth',clip(vScaledVec(inds(i))*5+0.25,0.5,inf));
        uistack(h,'bottom')
    end
end

% add a cone number to each cone
if withNum
    for i=1:size(dat.locs_c,1)
        h = text(dat.locs_c(i,1)+0.5,dat.locs_c(i,2)+0.5,num2str(dat.coneIds(i)));
        set(h,'Color','white');
    end
end

% wrangle the axes
set(gca,'YDir','reverse');
axis equal tight
xlim(dat.xRange);
ylim(dat.yRange);
axis off;