function plotModel2d(fits,celldat,dat,train,test,mode)

% plot the STA with dots

% plot the subunit fit
% - determine a centroid location
% - plot lines from subunit to centroid
% - draw an oval around each subunit

% - plot the output nonlinearity
% - plot the subunit nonlinearity
% - bar plot of r2

% plot the LN fit
% - use the same centroid location
% plot lines from cones to centroid

%% preamble stuff

if ~exist('newFig','var') || isempty(newFig)
    newFig = 1;
end

% set up a new window
if newFig
    clf
    set(gcf,'Position',[559 201 800 605]);
    set(gcf,'Color','w');
else
    cla
end

fit_SUB = fits.fit_SUB;
fit_LN = fits.fit_LN;
%fit_SING = fits.fit_SING;

out_SUB = evalFit(test,fit_SUB);
out_LN = evalFit(test,fit_LN);
%out_SING = evalFit(test,fit_SING);


%%
subselect = 1;
if subselect 
    diffVals = (out_SUB.Z_t - out_LN.Z_t).^2;
    inds = diffVals > prctile(diffVals,80) & diffVals < prctile(diffVals,100);
    %out_SING_SEL = evalFit(test,fit_SING,inds);
    out_SUB_SEL = evalFit(test,fit_SUB,inds);
    out_LN_SEL =  evalFit(test,fit_LN,inds);
    
    %diffVals2 = (out_SUB.Z_t - out_SING.Z_t).^2;
    %inds2 = diffVals2 > prctile(diffVals2,80) & diffVals2 < prctile(diffVals2,100);
    %out_SING_SEL2 = evalFit(test,fit_SING,inds2);
    %out_SUB_SEL2 = evalFit(test,fit_SUB,inds2);
    %out_LN_SEL2 =  evalFit(test,fit_LN,inds2);
end

%% plot the rf in pixel space

coneSz = 1.4;
subSz = 2;
coneEdge = 1.25;
figWidth = 0.3;
figHeight = 0.425;

clf

% plot the RF in pixel space
if isfield(celldat,'rf') && ~isempty(celldat.rf)
    h = axes('Position',[0.025 0.525 figWidth figHeight]);
    if celldat.polarity == -1
        image(1-norm_image(mean(celldat.rf,3)));
    elseif celldat.polarity == 1
        image(norm_image(mean(celldat.rf,3)));
    end
    xlim(celldat.xRange - celldat.xRange(1));
    ylim(celldat.yRange - celldat.yRange(1));
    axis image
    set(h,'YDir','reverse');
    axis off
end

% plot the cones for this cell
hold on
xPos = celldat.locs_c(:,1) - min(celldat.xRange) + 1;
yPos = celldat.locs_c(:,2) - min(celldat.yRange) + 1;
for i=1:length(xPos)
    circle(xPos(i),yPos(i),coneSz,'none','k',coneEdge);
end

% plot the other cones
if 1
hold on
[~,xInds] = isrange(dat.locations(:,1),...
  celldat.xRange(1),celldat.xRange(2));
[~,yInds] = isrange(dat.locations(:,2),...
  celldat.yRange(1),celldat.yRange(2));
inds = xInds & yInds;
inds(celldat.coneIds) = 0;
xPos = dat.locations(inds,1) - min(celldat.xRange) + 1;
yPos = dat.locations(inds,2) - min(celldat.yRange) + 1;
for i=1:length(xPos)
    circle(xPos(i),yPos(i),coneSz,'none','w',coneEdge);
end
end
%%

%%
h2 = axes('Position',[0.38 0.525 figWidth figHeight]);

axis equal
xlim(celldat.xRange);
ylim(celldat.yRange);
set(gca,'YDir','reverse');
axis off

rect1 = rectangle('Position',[celldat.xRange(1) celldat.yRange(1) ...
  sum(celldat.xRange) sum(celldat.yRange)]);
set(rect1,'FaceColor',[0.5 0.5 0.5]);
set(rect1,'EdgeColor','none');

%cell_border = drawConvexShape(dat.locs_c,[0.3 0.3 0.3],coneSz);
%set(cell_border,'EdgeColor','none');
plotSubunitsLinesOvals(fit_SUB,celldat,coneSz,subSz,coneEdge);

%%

%%

h3 = axes('Position',[0.75 0.6 0.2 0.2]);
plotNonLin(fit_SUB.g);
t = title('Output nonlinearity');
set(t,'FontSize',16);
set(t,'FontAngle','oblique');
set(gca,'FontSize',14);
%%

h3 = axes('Position',[0.75 0.775 0.2 0.2]);
plotNonLin(fit_SUB.f);
t = title('Subunit nonlinearity');
set(t,'FontSize',16);
set(t,'FontAngle','oblique');
set(gca,'FontSize',14);

%%
h3 = axes('Position',[0.38 0.05 figWidth figHeight]);
set(gca,'YDir','reverse');
axis equal
xlim(celldat.xRange);
ylim(celldat.yRange);
axis off
rect1 = rectangle('Position',[celldat.xRange(1) celldat.yRange(1)...
  sum(celldat.xRange) sum(celldat.yRange)]);
set(rect1,'FaceColor',[0.5 0.5 0.5]);
set(rect1,'EdgeColor','none');

centroid = getCentroid(celldat.locs_c);

% plot the lines
hold on
locs_s = normMat(fit_SUB.I_sc,'L1',2)*celldat.locs_c;
for ic=1:fit_SUB.C
    startPt = [celldat.locs_c(ic,1) celldat.locs_c(ic,2)];
    endPt = centroid;
    h = line([startPt(1) endPt(1)],[startPt(2) endPt(2)]);
    weight = (fit_LN.k(ic))/norm(fit_LN.k);   
    set(h,'Color','w');
    set(h,'LineWidth',abs(weight)*10);
end

% plot the cones
hold on
xPos = celldat.locs_c(:,1);
yPos = celldat.locs_c(:,2);
for i=1:size(celldat.locs_c,1)
    circle(xPos(i),yPos(i),coneSz,'w','k',coneEdge);
end

%%
h3 = axes('Position',[0.75 0.1 0.2 0.2]);
plotNonLin(fit_LN.g);
t = title('Output nonlinearity');
set(t,'FontSize',16);
set(t,'FontAngle','oblique');
set(gca,'FontSize',14);

%%

%axes('Position',[0.4 0.5 1 0.01]);
%h = line([0 1],[0.5 0.5]);
%set(h,'Color','k');
%axis off;

%%

h4 = axes('Position',[0.75 0.5 0.0825 0.1]);
R2 = [out_LN.r2; out_SUB.r2];
b = bar([1:2]-0.16,R2([1,2]),0.6);
set(b,'EdgeColor','none');
bar_child=get(b,'Children');
set(bar_child,'CData',[1 2]);
xlim([1-0.8 2+0.4]);
ylim([0 0.4]);
box off
set(gca,'YTick',[0 0.4]);
set(gca,'XTick',[]);
set(gca,'TickDir','out');
set(gca,'FontSize',14);
label = {'LN','Sub'};
for i=1:2
    t = text(i-0.14,-0.06,label{i});
    set(t,'HorizontalAlignment','center');
    set(t,'FontSize',14);
end
t = text(0.4,round2(max(R2),0.1,'ceil')+0.05,'All stimuli');
set(t,'FontSize',14,'FontAngle','italic');

h5 = axes('Position',[0.87 0.5 0.0825 0.1]);
R2_SEL = [out_LN_SEL.r2; out_SUB_SEL.r2];
b_SEL = bar([1:2]-0.14,R2_SEL([1,2]),0.6);
set(b_SEL,'EdgeColor','none');
bar_SEL_child=get(b_SEL,'Children');
set(bar_SEL_child,'CData',[1 2]);
uistack(b_SEL,'bottom');
colormap([0.75 0.75 0.75; 0.5 0.5 0.5; 0.25 0.25 0.25]);
xlim([1-0.8 2+0.4]);
ylim([0 0.4]);
box off
set(gca,'YTick',[]);
set(gca,'XTick',[]);
set(gca,'TickDir','out');
set(gca,'FontSize',14);
for i=1:2
    t = text(i-0.15,-0.06,label{i});
    set(t,'HorizontalAlignment','center');
    set(t,'FontSize',14);
end
t = text(0.6,round2(max(R2),0.1,'ceil')+0.05,'Subset');
set(t,'FontSize',14,'FontAngle','italic');


%%
v = loadVmat(celldat,dat);
if ~isempty(v)
    h5 = axes('Position',[0.025 0.05 figWidth figHeight]);
    rect1 = rectangle('Position',[celldat.xRange(1) celldat.yRange(1)...
      diff(celldat.xRange) diff(celldat.yRange)]);
    set(rect1,'FaceColor',[0.5 0.5 0.5]);
    set(rect1,'EdgeColor','none');
    plotPairwiseLines(celldat,v,0,0);
    uistack(rect1,'bottom');
    title(sprintf('Max linearity = %.02g',max(v(:))));
end

%%
plotTitles(celldat);
%%
set(gcf, 'InvertHardcopy', 'off');
set(gcf, 'PaperPositionMode', 'auto');
set(gcf,'PaperOrientation','landscape');
%saveName = sprintf('%g-%s-%s-2D',dat.rgcId,dat.dataSet(1:12),dat.dataSet(14:end));      
saveName = sprintf('%g-subunit-%s',celldat.rgcId,mode);
saveName = strcat(saveName,'-',num2str(celldat.percent));
savePath = sprintf('subunit-%s-%s-%s',mode,num2str(celldat.percent),celldat.cellType);

if ~isdir(fullfile(celldat.figurePath,savePath))
    mkdir(fullfile(celldat.figurePath,savePath))
end

print('-dpdf',fullfile(celldat.figurePath,savePath,strcat(saveName,'.pdf')));
