function plotPairwiseCell(dat,v,withNum)

%-------------------------------------------
%
% plotPairwiseCell(dat,v)
%
% plot aspects of a pairwise linearity measure 
% for a single ganlgion cell with C cones
%
% inputs:
% v -- a pairwise matrix of linearity, C x C
%      values on or below main diagonal should be 0
% 
% dat must contain the fields:
% dat.locs_c -- cone locations in xy, C x 2
% dat.xRange -- x range of the cone locations, 1 x 2
% dat.yRange -- y range of the cone locations, 1 x 2
%
% dat can contain the optional fields:
% dat.rf -- image of the sta
% dat.staFit -- weights to cones given by sta, C x 1
% dat.coneIds -- id number for each cone, C x 1
% dat.cellType -- string with cell type (e.g. 'off midget')
% dat.dataSet -- string with data set (e.g. 'date/piece/run')
% dat.rgcId -- rgc identifier number
%
% other optional inputs:
% withNum -- plot number next to each cone? (default = 0)
%
% freeman, 3-19-2012
%-------------------------------------------
%% preamble

if ~exist('withNum','var') || isempty(withNum)
    withNum = 0;
end

if withNum && ~isfield(dat,'coneIds');
    error('(plotPairwiseCell) must provide cone ids');
end

if isfield(dat,'rgcId') && isfield(dat,'cellType') && isfield(dat,'dataSet')
    titlePlot = 1;
else
    titlePlot = 0;
end

% set up the figure
figure
set(gcf,'Position',[800 80  600 730]);

%% plot the spatial sta if we have it
if isfield(dat,'rf');
    subplot(2,2,1);
    image(norm_image(dat.rf));
    axis equal
    axis off
end

%% plot the pairwise linearity in cone space
subplot(2,2,3);
colormap(cjet);
plotPairwiseLines(dat,v,withNum,1);

%% show the matrix of linearity
subplot(2,2,2);
imagesc(v); axis square off; colorbar;
if -max(v(:)) == max(v(:)) || -max(v(:)) > max(v(:))
set(gca,'CLim',[-0.1 0.1]);
else
set(gca,'CLim',[-max(v(:)) max(v(:))]);
end
title('Linearity index','FontSize',16);

%% plot the linearity measure as a function of distance
subplot(2,2,4);
D = getDistances(dat.locs_c,dat.locs_c);
vVec = triuVec(v,1);
dVec = triuVec(D,1);
if isfield(dat,'staFit');
    % get the average sta weight for each cone *pair*
    staMag = abs(dat.staFit)';
    staMag = staMag - min(staMag(:));
    staMag = staMag/max(staMag(:));
    staMagMat = (repmat(staMag,[size(staMag,2) 1]) + repmat(staMag,[size(staMag,2), 1])')/2;
    staMagMat = staMagMat - diag(diag(staMagMat));
    sVec = triuVec(staMagMat,1);
    [~,inds] = sort(sVec,'descend');
    if dat.polarity == -1
        h = scatter(dVec(inds),vVec(inds),30,-sVec(inds),'filled');
    else
        h = scatter(dVec(inds),vVec(inds),30,sVec(inds),'filled');
    end
    set(gca,'CLim',[-max(abs(sVec)) max(abs(sVec))]);
else
    h = scatter(dVec(inds),vVec(inds),30,'k','filled');
end
box off;
title('Linearity index versus distance','FontSize',16);
set(gca,'TickDir','out');
if isfield(dat,'staFit');
set(gca,'CLim',[-max(abs(sVec)) max(abs(sVec))]);
end
ylim([min(vVec(:))-2 max(vVec(:))+2]);

%% add titles to the plot
if titlePlot
    if withNum
        plotTitles(dat,'w');
    else
        plotTitles(dat,'k');
    end
end

%% set background to gray if we are using white on dark
if withNum
    set(gcf,'Color',[0.3 0.3 0.3]);
    set(gcf,'InvertHardcopy','off');
end