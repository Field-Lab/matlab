function save_2d_plots(xdata1, xdata2, y, mllparams_x, mllparams_y, cone1, cone2, sample)

if 0
    p = mllparams_x;
    
    [x1, x2] = meshgrid(-0.15:0.01:0.15,-0.15:0.01:0.15);
    sat   = p(1);
    sigma = p(2);
    mu = p(3);
    sh = p(4);
    a = p(5);
    b = p(6);
    y_interim = sat .* normcdf(a*x1 + b*x2, mu, sigma) + sh;
    
    figure;surface(x1, x2, y_interim);
    view(3)
    saveas(gcf, ['/Volumes/Analysis/2011-12-13-2/cone_data/manual/surface_plots_3736/cone', ...
        int2str(cone1), '_',int2str(cone2),'_sample_', int2str(sample),'_xshift_func.fig' ]);
    saveas(gcf, ['/Volumes/Analysis/2011-12-13-2/cone_data/manual/surface_plots_3736/cone', ...
        int2str(cone1), '_',int2str(cone2),'_sample_', int2str(sample),'_xshift_func.svg' ]);
    close all
    
    p = mllparams_y;
    sat   = p(1);
    sigma = p(2);
    mu = p(3);
    sh = p(4);
    a = p(5);
    b = p(6);
    y_interim = sat .* normcdf(a*x1, mu, sigma) + sat .* normcdf(b*x2, mu, sigma) + sh;
    
    figure;surface(x1, x2, y_interim);
    view(3)
    saveas(gcf, ['/Volumes/Analysis/2011-12-13-2/cone_data/manual/surface_plots_3736/cone', ...
        int2str(cone1), '_',int2str(cone2),'_sample_', int2str(sample),'_yshift_func.fig' ]);
    saveas(gcf, ['/Volumes/Analysis/2011-12-13-2/cone_data/manual/surface_plots_3736/cone', ...
        int2str(cone1), '_',int2str(cone2),'_sample_', int2str(sample),'_yshift_func.svg' ]);
    close all
end

clear fr
tmp = -0.15:0.02:0.15;
for i = 2:length(tmp)
    for j = 2:length(tmp)
        inds = find(xdata1>tmp(i-1) & xdata1<=tmp(i) & xdata2>tmp(j-1) & xdata2<=tmp(j));
        fr(i-1,j-1) = mean(y(inds));
    end
end
[x1_new, x2_new] = meshgrid(-0.14:0.02:0.15,-0.14:0.02:0.15);

figure;surface(x1_new, x2_new, fr);
view(3)
saveas(gcf, ['/Volumes/Analysis/2011-12-13-2/cone_data/manual/surface_plots_3736/data_only/data_', ...
    int2str(cone1), '_',int2str(cone2),'.svg' ]);
close all
