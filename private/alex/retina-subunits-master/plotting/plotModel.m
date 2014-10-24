function plotModel(fit,dat,nlPlot,titlePlot,newFig)

%-------------------------------------------
%
% plotModel(fit,dat,nlPlot,newFig)
%
% plot the components of an LN cascade model
% in diagram form
%
% fit must contain the fields:
% fit.S -- number of subunits
% fit.C -- number of stimuli
% fit.I_sc -- binary assignments from stimulus to subunits, S x C
% fit.A_sc -- weights from stimulus to subunits, S x C
% fit.J_s -- binary assignments from subunits to cell, 1 x S
% fit.B_s -- weights from subunits to cell, 1 x S

% optional fields for fit:
% fit.f -- parameters of subunit nonlinearity
% fit.g -- parameters of output nonlinearity
%
% optional input dat can contain the fields:
% dat.locs_c -- cone locations in xy, C x 2
% dat.rf -- image of the sta
% dat.cellType -- string with cell type (e.g. 'off midget')
% dat.dataSet -- string with data set (e.g. date/piece/run)
% dat.rgcId -- rgc identifier number
% dat.xRange -- x range of the cone locations, 1 x 2
% dat.yRange -- y range of the cone locations, 1 x 2
%
% more optional inputs:
% plotNonLin -- plot the nonlinearities? (default = 1);
% newFig -- make a new figure? (default = 1)
%
% freeman, 3-19-2012
%-------------------------------------------

%% preamble stuff
% check inputs
if ~exist('dat','var') || isempty(dat)
  dat = [];
end

if ~exist('nlPlot','var') || isempty(nlPlot)
  if isfield(fit,'f') && isfield(fit,'g')
    nlPlot = 1;
  else
    nlPlot = 0;
  end
end

if ~exist('titlePlot','var') || isempty(titlePlot)
  if isfield(dat,'rgcId') && isfield(dat,'cellType') && isfield(dat,'dataSet')
    titlePlot = 1;
  else
    titlePlot = 0;
  end
end

if ~exist('newFig','var') || isempty(newFig)
  newFig = 1;
end

% set up a new window
if newFig
  clf
  set(gcf,'Position',[800 80  600 730]);
else
  cla
end

% create some axes
axes('Position',[0.1 0.25 0.8 0.8]);

% specify some constants for sizes and locations
coneSz = 1;
subSz = 1.5;
rgcSz = 2;
conePosY = 0;
subPosY = -12.5;
rgcPosY = -20;
xPosMn = -10;
xPosMx = 10;
yBounds = [-24 5];
if nlPlot
  xBounds = [-11 20];
else
  xBounds = [-12 12];
end
widths = linspace(0.1,5,256);

%% plot the model diagram

% get colors for the subunits that have more than one cone
nonSingInds = sum(fit.I_sc,2) ~= 1;
nNonSing = sum(nonSingInds);
subClrs = zeros(fit.S,3);
subClrs(nonSingInds,:) = cmap_hsv_dark(nNonSing);

if isfield(dat,'locs_c') && ~isempty(dat.locs_c)
  % rescale to obtain a fixed ratio of the plot area
  % and the area containing the cones
  coneLocsX = dat.locs_c(:,1);
  coneLocsY = dat.locs_c(:,2);
  coneLocsX = (coneLocsX - mean(coneLocsX(:)));
  coneLocsY = (coneLocsY - mean(coneLocsY(:)));
  scaleFac2opts = linspace(0,5,100);
  for i=1:length(scaleFac2opts)
    tmpLocsX = coneLocsX*scaleFac2opts(i);
    tmpLocsY = coneLocsY*scaleFac2opts(i)/2;
    boundingBoxX = (max(tmpLocsX(:))-min(tmpLocsX(:)));
    boundingBoxY = (max(tmpLocsY(:))-min(tmpLocsY(:)));
    boundingArea = boundingBoxX*boundingBoxY;
    denRatio(i) = boundingArea/(xPosMx-xPosMn)^2;
  end
  [~,ind] = min(abs(denRatio-0.15));
  scaleFac2 = scaleFac2opts(ind);
  coneLocsX = coneLocsX * scaleFac2;
  coneLocsY = coneLocsY * scaleFac2/2;
  
else
  % if we don't have cone locations, just plot
  % rows of points
  coneLocsX = linspace(xPosMn,xPosMx,fit.C)-coneSz/2;
  coneLocsY = zeros(1,fit.C)+conePosY-coneSz/2;
end

% plot cones as circles
for ic=1:fit.C
  pos = [coneLocsX(ic) coneLocsY(ic)+coneSz/4 coneSz coneSz/2];
  rectangle('Position',pos,'Curvature',[1 1]);
end

% get subunit locations, sorted according to average x
% location of its cones (if locations are provided)
subLocsX = linspace(xPosMn,xPosMx,fit.S)-subSz/2;
subLocsY = zeros(1,fit.S)+subPosY-subSz/2;
if isfield(dat,'locs_c') && ~isempty(dat.locs_c)
  locs_s = normMat(fit.I_sc,'L1',2)*dat.locs_c;
  randn('seed',100);
  locs_s = locs_s + randn(size(locs_s))*0.0001;
  sorted_s = sort(locs_s(:,1),'ascend');
  for is=1:fit.S
    foo = find(sorted_s==locs_s(is,1));
    sub_inds_sorted(is) = foo(1);
  end
  subLocsX = subLocsX(sub_inds_sorted);
  subLocsY = subLocsY(sub_inds_sorted);
end

% scale the cone to subunit weights by dividing by max
A_sc_scaled = fit.A_sc;
A_sc_scaled = bsxfun(@rdivide,fit.A_sc',max(fit.A_sc'))';

% plot lines from cones to subunits, do it in order
% of the average y location of each subunit's cones (if provided),
% to give a nice 3d effect
if isfield(dat,'locs_c') && ~isempty(dat.locs_c)
  [~,sub_inds_sorted_y] = sort(locs_s(:,2),'ascend');
else
  sub_inds_sorted_y = [1:fit.S]';
end

for is=[sub_inds_sorted_y]'
  inds = find(fit.I_sc(is,:));
  for isc=1:length(inds)
    % get start and end points of lines from cone to subunit
    startPt = [coneLocsX(inds(isc)) coneLocsY(inds(isc))]+coneSz/2;
    endPt = [subLocsX(is) subLocsY(is)]+subSz/2;
    % plot the line, with weight/color depending on
    % whether it's a singleton or not
    h = line([startPt(1) endPt(1)],[startPt(2) endPt(2)]);
    weight = ceil(A_sc_scaled(is,inds(isc))*254)+1;
    if nonSingInds(is) == 1
      set(h,'LineWidth',widths(weight));
      set(h,'Color',subClrs(is,:));
      uistack(h,'top');
    else
      set(h,'LineWidth',0.5);
      set(h,'Color','k');
      uistack(h,'bottom');
    end
  end
end

% plot the subunits as circles
for is=1:fit.S
  pos = [subLocsX(is) subLocsY(is)+subSz/4 subSz subSz/2];
  if subClrs(is,:) == [0 0 0];
    tmpClr = [1 1 1];
  else
    tmpClr = subClrs(is,:);
  end
  rectangle('Position',pos,'Curvature',[1 1],'FaceColor',tmpClr);
end

% plot weighted lines from subunits to rgc
rgcLocX = mean([xPosMx xPosMn]) - rgcSz/2;
rgcLocY = rgcPosY - rgcSz/2;
B_s_scaled = fit.B_s;
B_s_scaled = B_s_scaled/max(abs(B_s_scaled));
B_s_scaled_mag = abs(B_s_scaled);
for is=1:fit.S
  startPt = [subLocsX(is) subLocsY(is)]+subSz/2;
  endPt = [rgcLocX(1) rgcLocY(1)]+rgcSz/2;
  h = line([startPt(1) endPt(1)],[startPt(2) endPt(2)]);
  weightClr = clip(round(B_s_scaled(is)*256/2+256/2),1, 256);
  weightLine = clip(round(B_s_scaled_mag(is)*256),1,256);
  set(h,'LineWidth',widths(weightLine));
  set(h,'Color','k');
end

% plot the rgc as a circle
pos = [rgcLocX(1) rgcLocY(1) rgcSz rgcSz];
rectangle('Position',pos,'Curvature',[1 1],'FaceColor',[1 1 1]);

% do some axis wrangling
axis equal
xlim(xBounds);
ylim(yBounds);
axis off
set(gcf,'Color','white')
set(gca,'YDir','reverse');

%% plot the nonlinearities
if nlPlot
  axes('Position',[0.75 0.58 0.2 0.2]);
  plotNonLin(fit.f);
  title('Subunit nonlinearity');
  axes('Position',[0.75 0.775 0.2 0.2]);
  plotNonLin(fit.g);
  title('Output nonlinearity');
end

%% plot the rf in pixel space
if isfield(dat,'rf') && ~isempty(dat.rf)
  axes('Position',[0.75 0.35 0.2 0.2]);
  image(norm_image(dat.rf));
  axis off equal
  set(gca,'YDir','reverse');
end

%% add titles to the plot
if titlePlot
  plotTitles(dat);
end

%% plot each subunit separately in cone space
if isfield(dat,'locs_c') && ~isempty(dat.locs_c)
  
  % figure out how many plots we need, and arrange them
  % to tile as squares inside a rectangle
  inds = find(sum(fit.I_sc,2)>1);
  nPlots = length(inds);
  left = 0.05;
  bottom = 0.075;
  horzWidth = 0.9;
  vertWidth = 0.25;
  
  % do we need 1 row of plots or 2?
  if horzWidth/nPlots < vertWidth/2; nRows = 2; else nRows = 1; end
  if nRows == 1
    sz = min(vertWidth,horzWidth/nPlots);
    xPos = left + [0:nPlots-1]*sz;
    yPos = ones(1,nPlots)*(bottom);
  end
  if nRows == 2
    sz = min(vertWidth/2,horzWidth/ceil(nPlots/2));
    nPlotsWidth = floor(horzWidth/sz);
    xPos = left + [0:nPlotsWidth-1]*sz;
    xPos = [xPos, xPos];
    yPos1 = ones(1,nPlotsWidth)*(bottom+sz);
    yPos2 = ones(1,nPlotsWidth)*bottom;
    yPos = [yPos1, yPos2];
  end
  % loop over subunits
  for is=1:length(inds)
    % create an axes
    a = axes('Position',[xPos(is) yPos(is) sz sz]);
    % get the weights from cones to this subunit
    weights = fit.A_sc(inds(is),:);
    % use the same color for each subunit as from above,
    % but adjust the saturation of each cone's
    % color according to the weight
    clrsThisSub = repmat(subClrs(inds(is),:),[length(weights) 1]);
    clrsThisSub_hsv = rgb2hsv(clrsThisSub);
    clrsThisSub_hsv(:,2) = clrsThisSub_hsv(:,2).*weights(:)/max(weights(:));
    clrsThisSub_hsv(:,3) = clrsThisSub_hsv(:,3)+0.2*(1-weights(:)/max(weights(:)));
    clrsThisSub = hsv2rgb(clrsThisSub_hsv);
    clrsThisSub(weights(:)==0,1) = 1;
    clrsThisSub(weights(:)==0,2) = 1;
    clrsThisSub(weights(:)==0,3) = 1;
    % plot the individual cones
    % (use this instead of scatter because scatter
    % fucks stuff up in the rest of the plot)
    hold on
    for ic=1:length(weights)
      h = plot(dat.locs_c(ic,1),dat.locs_c(ic,2),'o');
      set(h,'MarkerFaceColor',clrsThisSub(ic,:));
      set(h,'MarkerSize',8);
      set(h,'MarkerEdgeColor','k');
    end
    hold off
    % wrangle the axes
    xlim(dat.xRange+[-2 2]);
    ylim(dat.yRange+[-2 2]);
    uistack(a,'bottom');
    axis equal off;
    set(gca,'CLim',[-max(weights(:)) max(weights(:))]);
    set(gca,'YDir','reverse');
  end
end