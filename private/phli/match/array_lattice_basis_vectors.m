%%
ai = load_array_info(1504);
positions = [ai.positions'; ones(1, size(ai.positions, 1))];
basis = [30 0 0; 15 30 0; 0 0 1];
translation = [1 0 -13; 0 1 -14; 0 0 1];
T = (basis*translation) \ positions;
TT = int16(T(1:2,:))';


%% Test; fails for some because load_array_info is returning Xe-15 instead of 0!
basis*translation*T == positions


%%
subplot(1, 2, 1);
plot(TT(:,1), TT(:,2), '.')
for i = 1:size(TT, 1)
    text(double(TT(i,1)), double(TT(i,2)), num2str(i));
end

subplot(1, 2, 2);
plot(ai.positions(:,1), ai.positions(:,2), '.')
for i = 1:size(ai.positions, 1)
    text(ai.positions(i,1), ai.positions(i,2), num2str(i));
end
clear i


%% Simple latticed test
L = nan(max(TT));
ind = sub2ind(size(L), TT(:,1), TT(:,2));
L(ind) = 1;
clear ind;
surf(L);


%% Load EI
datarun = load_data('2011-07-14-6/data005');
datarun = load_params(datarun);
datarun = load_ei(datarun, []);
ei = get_ei_max_frame(get_ei(datarun, datarun.cell_types{1}.cell_ids(4)));


%% Latticed EI test
L = nan(max(TT));
ind = sub2ind(size(L), TT(:,1), TT(:,2));
L(ind) = abs(ei);
clear ind;
surf(L);


%%
p0 = 37/20; p1 = -41/240; p2 = 7/240;
prefilt = [
    0  p2 0  0  0;
    p2 p1 p1 p2 0;
    0  p1 p0 p1 0;
    0  p2 p1 p1 p2;
    0  0  0  p2 0;
];
clear p0 p1 p2;


%% Get 3-dir box spline quasi-interpolation coefficients
L = zeros(max(TT));
ind = sub2ind(size(L), TT(:,1), TT(:,2));
L(ind) = abs(ei);
clear ind;
C = conv2(prefilt, L);
surf(C);


%% Get coarse 2nd order boxspline
xi = 1:9;
yi = xi;
[XI, YI] = meshgrid(xi, yi);
error('I think this should be NDGRID instead of MESHGRID!');

bstranslation = [-5 -5];
hexbasis = [1 0.5; 0 sqrt(3)/2];
hexpoints = hexbasis*[XI(:)+bstranslation(1) YI(:)+bstranslation(2)]';

bs2vals = boxSplineDn(hexpoints(1,:), hexpoints(2,:), 2);
B2 = zeros(size(XI));
ind = sub2ind(size(XI), XI(:), YI(:));
B2(ind) = bs2vals;
clear xi yi ind;


%% Reconstruct L
LL = conv2(C, B2);
surf(LL);