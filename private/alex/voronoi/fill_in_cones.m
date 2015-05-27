datarun = load_data('/Volumes/Analysis/2011-12-13-2/d08-11-norefit/data008-from-d08_11/data008-from-d08_11');
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun);
datarun = load_cones_ath(datarun, 'd08-11-norefit_d08-bayes-msf_10.00');


x = datarun.cones.centers(:,1);
y = datarun.cones.centers(:,2);
ncones = length(x);

figure
plot(x, y, '*')
axis ij

dists = squareform(pdist([x y]));
dists(dists==0) = NaN;

tmp = sort(dists);
tmp = tmp(1:3,:);
figure
hist(tmp(:), 1:1:20)
tmp = tmp(tmp<15);
mean_nnd = robust_mean(tmp(:));
mean_nnd_std = robust_std(tmp(:));

tmp = sort(dists);
tmp = tmp(1:6,:);

for i = 1:ncones
    [tmp, ic] = sort(dists(:,i));
    close_cones = find(tmp(1:6)<mean_nnd+2*mean_nnd_std);
    neighbours = ic(close_cones);
    squareform(pdist([x(neighbours)  y(neighbours)]))

end



[val, inds] = min(dists);
bins = [0:.5:10 20];
figure
hist(val,bins);

mean_space = robust_mean(val);
mean_space_std = robust_std(val);



figure
plot(x, y, '*')
axis ij

full_map = zeros(300,300);
linind = sub2ind([300,300],round(y), round(x));
full_map(linind) = 1:ncones;
figure
imagesc(full_map)
set(gca,'dataaspectratio',[1 1 1])

% find holes in mosaic?
% first go through existing cones. Check out each pair of closest neighbours.

for i = 1:ncones
    if dists(i,inds(i))<(mean_space+3*mean_space_std)
        

    end
end


%%
weaker = 5;
% func: f(x) = a/x; a = 6.3
mean_space = 6.2;
sp = 6.3/mean_space;
bords = round(mean_space)*3;

d = (-bords:1/3:bords);  % domain for both x,y
f = @(x,y) (cos(sqrt((x*sp).^2 + (y*sp).^2)).*normpdf(x,0,weaker).*normpdf(y,0,weaker)); % a function of two variables                   
[X,Y] = meshgrid(d,d);                  % create a grid of points
Z = f(X,Y);                             % evaluate f at every point on grid
Z = Z/max(abs(Z(:)));

plot(d,Z(:,round(end/2)))
imagesc(Z)

figure
plot((76-[6 13 24 26 28 41 45 51 55 60])/10,[0.9 1 1.2 1.26 1.3 1.8 2 2.5 3 4])


a = (76-[6 13 24 28 41 45 51 55 60])/10;
b = [0.9 1 1.2 1.3 1.8 2 2.5 3 4];



%%
x = round(datarun.cones.centers(:,1)*3);
y = round(datarun.cones.centers(:,2)*3);


c1 = 50;
c2 = 46;
c3 = 47;

full_map = zeros(300*3);
full_map2 = zeros(300*3);
full_map3 = zeros(300*3);
full_map(x(c1)-bords*3:x(c1)+bords*3, y(c1)-bords*3:y(c1)+bords*3) = Z;
full_map2(x(c2)-bords*3:x(c2)+bords*3, y(c2)-bords*3:y(c2)+bords*3) = Z;
full_map3(x(c3)-bords*3:x(c3)+bords*3, y(c3)-bords*3:y(c3)+bords*3) = Z;
figure
imagesc(full_map)
figure
imagesc(full_map2)
figure
imagesc(full_map3)

figure
imagesc(full_map+full_map2)%+full_map3);

%% Fourier
img   = full_map
imagesc(img)
F     = fft2(img);

figure;

imagesc(100*log(1+abs(fftshift(F)))); colormap(gray); 
title('magnitude spectrum');

figure;
imagesc(angle(F));  colormap(gray);
title('phase spectrum');




