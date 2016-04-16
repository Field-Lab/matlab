function save_2d_plots(xdata1, xdata2, y, mllparams_x, mllparams_y, cone1, cone2, sample)

if 1
    p = mllparams_x;
    
    tt = min(max(xdata1), max(xdata2));
    tt1 = max(min(xdata1), min(xdata2));
    tmp = linspace(tt1, tt, 7);
    
    [x1, x2] = meshgrid(tmp,tmp);
    
    sat   = p(1);
    sigma = p(2);
    mu = p(3);
    sh = p(4);
    a = p(5);
    y_interim = sat .* normcdf(a*x1 + x2, mu, sigma) + sh;
    
    figure;surface(x1, x2, y_interim);
    view(3)
    title('x shift')

    
    p = mllparams_y;
    sat   = p(1);
    sigma = p(2);
    mu = p(3);
    sh = p(4);
    a = p(5);
    y_interim = sat .* (normcdf(a*x1, mu, sigma) + normcdf(x2, mu, sigma)) + sh;
    
    figure;surface(x1, x2, y_interim);
    view(3)
    title('y shift')
end

clear fr

tt = min(max(xdata1), max(xdata2));
tt1 = max(min(xdata1), min(xdata2));
tmp1 = linspace(tt1, tt, 7);
tmp2 = tmp1;
% tmp1 = linspace(min(xdata1), max(xdata1), 10);
% tmp2 = linspace(min(xdata2), max(xdata2), 10);
for i = 2:length(tmp1)
    for j = 2:length(tmp2)
        inds = find(xdata1>tmp1(i-1) & xdata1<=tmp1(i) & xdata2>tmp2(j-1) & xdata2<=tmp2(j));
        fr(i-1,j-1) = mean(y(inds));
    end
end
[x1_new, x2_new] = meshgrid(tmp1(1:end-1),tmp2(1:end-1));

figure;surface(x1_new, x2_new, fr);
view(3)
% saveas(gcf, ['/Volumes/Analysis/2016-02-17-4/cone_data/manual/surface_plots_3274/data_only/data_', ...
%     int2str(cone1), '_',int2str(cone2),'.bmp' ]);
% close all

% saveas(gcf, ['/Volumes/Analysis/2011-12-13-2/cone_data/manual/surface_plots_3736/data_only/data_', ...
%     int2str(cone1), '_',int2str(cone2),'.svg' ]);
% close all
