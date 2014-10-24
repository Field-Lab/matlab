function plotSTC(dat,fit,k)

%-------------------------------------------
%
% plotSTC(data,fit)
%
% plot the results of an STC analysis
% fit should contain 
%
% inputs:
% data -- stucture containing the data
% data.locs_s -- locations of cones
%
% fit -- structure containing fit parameters
% fit.k -- sta
% fit.V -- matrix of stc filters
% fit.e -- eigenvalues
% fit.e_rnd -- eigenvalues after randomization
% k -- number of filters to plot
%
% freeman, 4-1-2012
%-------------------------------------------

if ~exist('k','var') || isempty(k)
    k = 6;
end

if k > 6
    fprintf('(plotSTC) can only plot up to 6 filters');
    k = 6;
end

if size(fit.V,2) < k
    k = size(fit.V,2);
end

figure;
set(gcf,'Position',[800 80  600 730]);
colormap(cjet);

if isfield(dat,'rf') && ~isempty(dat.rf)
subplot(3,3,1);
image(norm_image(dat.rf));
axis equal
axis off
end

subplot(3,3,2);
scatter(dat.locs_c(:,1),dat.locs_c(:,2),50,fit.k(:),'filled');
set(gca,'CLim',[-max(abs(fit.k(:))) max(abs(fit.k(:)))]);
set(gca,'YDir','reverse');
xlim(dat.xRange);
ylim(dat.yRange);
axis equal
axis off;
title('STA','FontSize',14);

subplot(3,3,3);
mxEigN = 20;
if length(fit.e) < mxEigN; mxEigN = length(fit.e); end;
if isfield(fit,'e_rnd');
    hold on
    plot(prctile(fit.e_rnd(:,1:mxEigN),[1]),'r','LineWidth',3);
    plot(prctile(fit.e_rnd(:,1:mxEigN),[99]),'r','LineWidth',3);
end
h = plot(fit.e(1:mxEigN),'k.');
set(h,'MarkerSize',20);
xlim([0 mxEigN+1]);
set(gca,'TickDir','out');
title('Eigenvalue spectrum','FontSize',14);
box off;

ic = 4;
for ik=1:k
    subplot(3,3,ic);
    scatter(dat.locs_c(:,1),dat.locs_c(:,2),50,fit.V(:,ik),'filled');
    set(gca,'CLim',[-max(abs(vector(fit.V(:,ik)))) max(abs(vector(fit.V(:,ik))))]);
    set(gca,'YDir','reverse');
    xlim(dat.xRange);
    ylim(dat.yRange);
    ic = ic + 1;
    axis off;
    axis equal;
    title(sprintf('STC %g',ik),'FontSize',14);
end

plotTitles(dat,'k');
