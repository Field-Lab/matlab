
% to create a voronoi map with rings:
% 1. Identify cone centers (see above).
% 2. Plot a certain shape (1, cross-like 5, or square 3x3 pixels) around
% the center. this is the main voronoi region.
% 3. Plot a ring around the main Voronoi region, with given width (1,2
% pixels... depends on cone spacing)
% 4. Fill in remaining space with single pixel noise (or more. Depends
% on cone spacing).

% 5. Overlaps: look for closest cone center (floating point!) and make it
% the same. If equal, choose randomly.

% check for missing cones! ask Nishal


datarun = load_data('/Volumes/Acquisition/Analysis/2015-10-06-5/data006/data006');
datarun = load_params(datarun,'verbose',1);
% datarun = load_cones_ath(datarun, 'd08-11-norefit_d08-bayes-msf_10.00');

% datarun = load_sta(datarun);
% datarun = set_polarities(datarun);
% datarun = load_neurons(datarun);

datarun.cones.centers = stim.coord;
figure
plot(datarun.cones.centers(:,1), datarun.cones.centers(:,2), '*')
axis ij
x = round(datarun.cones.centers(:,1));
y = round(datarun.cones.centers(:,2));
linind = sub2ind([320,320],y, x);
ncones = size(datarun.cones.centers,1);



sidesize = 3; % side of the square, or diameter of the circle. Should be odd.
ringsize = 0;
shape = 'round'; % options are 'square','cross', 'round'


full_map = zeros(320);
full_map(linind) = 1:ncones;

% figure the shifts of linear indices
switch shape
    case 'square'
        a = ones(sidesize);
    case 'cross'
        a = ones(sidesize);
        a(1,1) = 0;
        a(1,end) = 0;
        a(end,1) = 0;
        a(end,end) = 0;
    case 'round'
        [rr cc] = meshgrid(1:sidesize+1);
        rad = (sidesize+1)/2;
        C = sqrt((rr-rad).^2+(cc-rad).^2)<=sidesize/2;
        a = double(C(1:end-1,1:end-1)); 

end

main_map = conv2(full_map,a,'same');

figure
imagesc(main_map)
set(gca, 'dataaspectratio',[1 1 1])
hold on
plot(x,y,'.r')

% tmp = imresize(main_map,0.25,'nearest');
% figure
% imagesc(tmp)
% set(gca, 'dataaspectratio',[1 1 1])

% make rings (as a separate map) - map 1;
full_map = zeros(320);
full_map(linind) = 1:ncones;
[rr cc] = meshgrid(1:ringsize+1);
rad = (ringsize+1)/2;
C = sqrt((rr-rad).^2+(cc-rad).^2)<=ringsize/2;
a = double(C(1:end-1,1:end-1));
ring_map = conv2(full_map,a,'same');

% make rings (as a separate map) - map 2;

full_map = zeros(320);
full_map(linind) = 2:ncones+1;
[rr cc] = meshgrid(1:ringsize+1);
rad = (ringsize+1)/2;
C = sqrt((rr-rad).^2+(cc-rad).^2)<=ringsize/2;
a = double(C(1:end-1,1:end-1));
ring_map2 = conv2(full_map,a,'same');


% compare 2 maps
m1 = ring_map;
m2 = ring_map2;

m2(m2>0) = m2(m2>0)-1;
a = m1-m2;
m3 = m1;
m3(find(a))=ncones+1;

tmp = m3;
tmp(tmp>0) = tmp(tmp>0)+ncones;
tmp(main_map>0) = 0;
final_map = main_map+tmp;

figure
imagesc(final_map)
set(gca, 'dataaspectratio',[1 1 1])

linind = find(final_map==max(final_map(:)));
[a,b] = find(final_map==max(final_map(:)));

tmp = pdist2([a, b], [y, x]);
[val, ind] = min(tmp');

final_map(linind)=ind + ncones;
figure
imagesc(final_map)
set(gca, 'dataaspectratio',[1 1 1])


% segments

tmp_map = final_map;

for i = 1:length(x)
    my_val = final_map(y(i),x(i))+ncones;
    [a, b] = find(final_map==my_val);
    lins = sub2ind([320,320],a, b);
    
    my_nums = a<=y(i)&b<=x(i);
    tmp_map(lins(my_nums)) = my_val;
    
    my_nums = a>=y(i)&b<=x(i);
    tmp_map(lins(my_nums)) = my_val+ncones;
    
    my_nums = a<=y(i)&b>=x(i);
    tmp_map(lins(my_nums)) = my_val+ncones*2;
    
    my_nums = a>=y(i)&b>=x(i);
    tmp_map(lins(my_nums)) = my_val+ncones*3;
  
end

final_map = tmp_map;
figure
imagesc(final_map)
set(gca, 'dataaspectratio',[1 1 1])


% final_map = main_map; % if rings are 0
%% save stuff
my_info.datarun = '/Volumes/Analysis/2015-10-06-2/data004/data004';
my_info.cone_map = 'manual_s_cones_map004_NOrings_small';
my_info.sidesize = sidesize; 
my_info.ringsize = ringsize;
my_info.shape = shape;
my_info.segments = 4; % set the number of segments! (currently 1 or 4)

my_dir = '/Volumes/Analysis/2015-10-06-2/maps/data004_no_rings_small/';

mkdir(my_dir);
dlmwrite([my_dir 'map-NOrings-small-004.txt'], final_map, 'delimiter', '\t', 'newline', 'pc');
save([my_dir '/info'],'my_info')

% check
final_map = load('/Volumes/Analysis/2015-10-06-2/maps/data004_rings/map-rings-004.txt');
figure; imagesc(final_map);
load('/Volumes/Analysis/2015-10-06-2/maps/info.mat')
