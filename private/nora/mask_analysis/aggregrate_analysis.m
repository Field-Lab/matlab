default_colors = get(gca,'ColorOrder');
thresh = 0;
exp = 1;
load('SpotStats2016-01-05-01.mat')
NSEM_Corr = NSEM_Corr./repmat(max(NSEM_Corr'), [2,1])';
NSEM_Corr(sum(diff(NSEM_Corr')>thresh)>0,:) = 0;
on = [1 3];
off = [2 4];
sigmas = [2 4];
hold on;
plot(sigmas, NSEM_Corr(on, :)', 'Color', default_colors(exp,:), 'LineWidth', 2);
plot(sigmas, NSEM_Corr(off, :)', 'Color', default_colors(exp,:));

load('SpotStats2016-01-05-02.mat')
NSEM_Corr = NSEM_Corr./repmat(max(NSEM_Corr'), [2,1])';
NSEM_Corr(sum(diff(NSEM_Corr')>thresh)>0,:) = 0;
on = [1 3];
off = [2 4];
sigmas = [2 4];
plot(sigmas, NSEM_Corr(on, :)', 'Color', default_colors(exp,:), 'LineWidth', 2);
plot(sigmas, NSEM_Corr(off, :)', 'Color', default_colors(exp,:));

exp = 2;
load('SpotStats2016-02-17-1data021.mat')
NSEM_Corr = NSEM_Corr./repmat(max(NSEM_Corr'), [4,1])';
NSEM_Corr(sum(diff(NSEM_Corr')>thresh)>0,:) = 0;
on = 1:3;
off = 4:6;
sigmas = [2 4 5 6];
plot(sigmas, NSEM_Corr(on, :)','Color', default_colors(exp,:), 'LineWidth', 2);
plot(sigmas, NSEM_Corr(off, :)', 'Color', default_colors(exp,:));

% load('SpotStats2016-02-17-1data023.mat')
% NSEM_Corr(sum(diff(NSEM_Corr')>-0.2)>0,:) = 0;
% on = 1:2:6;
% off = 2:2:6;
% sigmas = [2 4 5 6];
% plot(sigmas, NSEM_Corr(on, :)','--', 'Color', default_colors(exp,:), 'LineWidth', 2);
% plot(sigmas, NSEM_Corr(off, :)', 'Color', default_colors(exp,:));

exp = 3;
load('SpotStats2016-02-17-6data002.mat')
NSEM_Corr = NSEM_Corr./repmat(max(NSEM_Corr'), [4,1])';
NSEM_Corr(sum(diff(NSEM_Corr')>thresh)>0,:) = 0;
off = 1:3;
on = 4:6;
sigmas = [2 4 5 6];
plot(sigmas, NSEM_Corr(on, :)','Color', default_colors(exp,:), 'LineWidth', 2);
plot(sigmas, NSEM_Corr(off, :)', 'Color', default_colors(exp,:));

load('SpotStats2016-02-17-6data004.mat')
NSEM_Corr = NSEM_Corr./repmat(max(NSEM_Corr'), [4,1])';
NSEM_Corr(sum(diff(NSEM_Corr')>thresh)>0,:) = 0;
on = 1:3; 
off = 4:6;
sigmas = [2 4 5 6];
plot(sigmas, NSEM_Corr(on, :)','Color', default_colors(exp,:), 'LineWidth', 2);
plot(sigmas, NSEM_Corr(off, :)', 'Color', default_colors(exp,:));

exp = 4;
load('SpotStats2016-02-17-8data002.mat')
NSEM_Corr = NSEM_Corr./repmat(max(NSEM_Corr'), [4,1])';
NSEM_Corr(sum(diff(NSEM_Corr')>thresh)>0,:) = 0;
on = 1:2:5;
off = 2:2:6;
sigmas = [2 4 5 6];
plot(sigmas, NSEM_Corr(on, :)','Color', default_colors(exp,:), 'LineWidth', 2);
plot(sigmas, NSEM_Corr(off, :)', 'Color', default_colors(exp,:));
