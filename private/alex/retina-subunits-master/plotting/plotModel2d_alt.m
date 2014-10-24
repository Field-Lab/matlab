function plotModel2d_alt(fits,celldat,dat,train,test,mode)

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
    set(gcf,'Position',[559 201 800 303]);
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
end

%% plot the rf in pixel space

coneSz = 1.4;
subSz = 2;
coneEdge = 1.25;
figWidth = 0.3;
figHeight = 0.85;

clf

% plot the RF in pixel space
if isfield(celldat,'rf') && ~isempty(celldat.rf)
    h = axes('Position',[0.025 0.07 figWidth figHeight]);
    image(norm_image(mean(celldat.rf,3)));
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


%%
h2 = axes('Position',[0.37 0.07 figWidth figHeight]);

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

h3 = axes('Position',[0.68 0.51 0.35 0.35]);
plotNonLin(fit_SUB.f);
t = title('Subunit nonlinearity');
set(t,'FontSize',16);
set(t,'FontAngle','oblique');
h = xlabel('Input');
set(h,'FontSize',14);
set(h,'VerticalAlignment','baseline');
set(h,'FontAngle','italic');
h = ylabel('Output');
set(h,'FontSize',14);
set(h,'VerticalAlignment','baseline');
set(h,'FontAngle','italic');
set(gca,'FontSize',14);


h4 = axes('Position',[0.75 0.1 0.21 0.2]);
R2 = [out_LN.r2; out_SUB.r2];
if subselect
    b = bar([1:length(R2)]-0.16,R2,0.25);   
else
    b = bar([1:length(R2)],R2,0.25);
end
set(b,'EdgeColor','none');
bar_child=get(b,'Children');
set(bar_child,'CData',[1 2]);
if subselect
    hold on
    R2_SEL = [out_LN_SEL.r2; out_SUB_SEL.r2];
    b_SEL = bar([1:length(R2_SEL)]+0.16,R2_SEL,0.25);
    set(b_SEL,'EdgeColor','none');
    bar_SEL_child=get(b_SEL,'Children');
    set(bar_SEL_child,'CData',[1 3]);
    uistack(b_SEL,'bottom');
    hold off
end

colormap([0.75 0.75 0.75; 0.25 0.25 0.25]);
ylim([0 0.5]);
xlim([1-0.5 length(R2)+0.5]);
box off
t = title('Model performance');
set(t,'FontSize',16);
set(t,'FontAngle','oblique');
set(gca,'YTick',[0 0.5]);
set(gca,'XTick',[]);
set(gca,'FontSize',14);
label = {'LN','Subunit'};
for i=1:2
    t = text(i,-0.1,label{i});
    set(t,'HorizontalAlignment','center');
    set(t,'FontSize',14);
end


set(gca,'TickDir','out');

for i=2:2
    %if subselect
    %    t = text(i-0.14,R2(i)+0.06,strcat(sprintf('%.02f',100*(R2(i)/R2(1)-1)),'%'));
    %else
    %    t = text(i,R2(i)+0.06,strcat(sprintf('%.02f',100*(R2(i)/R2(1)-1)),'%'));
    %end
    set(t,'HorizontalAlignment','center');
    set(t,'FontSize',14);
    set(t,'FontAngle','italic');
end
t = text(0.1, 0.25, 'R2');
set(t,'FontSize',14);


%%
plotTitles(celldat);