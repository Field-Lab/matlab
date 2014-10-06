%% Calculate rod spacing from density and hexagonal lattice
rods_per_1deg_diameter = 10000;
deg2_per_rod = 0.5*0.5*pi/rods_per_1deg_diameter;
rod_hex_side_length = 2*sqrt(2/3/sqrt(3)*deg2_per_rod);
rod_spacing = 2 * rod_hex_side_length;


%% Build pools of 4-cycle rods

% Get coordinates of seed rods
xi = -20:2:20;
yi = xi;
[XI,YI] = ndgrid(xi,yi);

hex_basis = rod_spacing .* [
    0  sqrt(3)/2;
    1  1/2;
];

seed_positions = hex_basis * [XI(:) YI(:)]';

% Build pools off of seed positions
for i = 1:size(seed_positions, 2)
    spx = seed_positions(1,i);
    spy = seed_positions(2,i);
    poolx(i,:) = [spx+hex_basis(1,1) spx spx+hex_basis(1,2) spx+hex_basis(1,1)+hex_basis(1,2)];
    pooly(i,:) = [spy+hex_basis(2,1) spy spy+hex_basis(2,2) spy+hex_basis(2,1)+hex_basis(2,2)];
    
    % For plotting, nice to add an extra point to close the cycle
%     poolx(i,:) = [spx+hex_basis(1,1) spx spx+hex_basis(1,2) spx+hex_basis(1,1)+hex_basis(1,2) spx+hex_basis(1,1)];
%     pooly(i,:) = [spy+hex_basis(2,1) spy spy+hex_basis(2,2) spy+hex_basis(2,1)+hex_basis(2,2) spy+hex_basis(2,1)];
end


%% Get normalized illumination of pools

% 2.5mm pupil, Artal & Navarro 1993
% A = 0.16; B = 0.06; C = 0.36;
% savename = 'pool_illumination_2_5mm';

% 8mm pupil, Artal & Navarro 1993
% A = 0.53; B = 0.08; C = 0.11;

% 3mm pupil, Guirao et al. 1999
a = 16.12; b = 17.5;
A = 1/a; B = 1/b; C = 1/4;
savename = 'pool_illumination_3mmb';

pool_illumination = zeros(size(poolx));
for i = 1:numel(poolx)
    x = poolx(i); y = pooly(i);
    d = sqrt(x^2 + y^2);
    pool_illumination(i) = (1-C)*2*A./(A*A + 4*pi*pi*d*d) + C*2*B./(B*B + 4*pi*pi*d*d);
end

total_illumination = sum(sum(pool_illumination));
pool_illumination = pool_illumination ./ total_illumination;

plot3(poolx', pooly', pool_illumination'); set(gca, 'DataAspectRatio', [1 1 0.02]);

save(savename, 'poolx', 'pooly', 'pool_illumination');

%% Looks right...
d = sqrt(poolx.^2 + pooly.^2);
plot(d, total_illumination .* pool_illumination, 'ko');