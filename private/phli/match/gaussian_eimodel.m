function eimodel = gaussian_eimodel(params, times, positions)

mux    = params(1);
muy    = params(2);
a      = params(3);
sigma  = params(4);
mut    = params(5);
sigmat = params(6);


dist = sqrt((positions(:,1)-mux).^2 + (positions(:,2)-muy).^2);
space = a * normpdf(dist, 0, sigma);
time = normpdf(times(:)', mut, sigmat);

eimodel = space * time;