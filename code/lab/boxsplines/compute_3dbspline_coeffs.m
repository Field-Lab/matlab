function [cx, cy, coeffs] = compute_3dbspline_coeffs(x, y, z, bases, prefilt)
% Default prefilt is for 2nd order 3dbsplines.

if nargin < 4
    bases = [30 0; 15 30];
end

if nargin < 5
    p0 = 37/20; p1 = -41/240; p2 = 7/240;
    prefilt = [
        0  p2 0  0  0;
        p2 p1 p1 p2 0;
        0  p1 p0 p1 0;
        0  p2 p1 p1 p2;
        0  0  0  p2 0;
    ];
    clear p0 p1 p2;
end


% Calculate coordinates in new basis
rebased_coords = bases \ [x(:)'; y(:)'];


% Shift to get only positive indices
trans = abs(min(rebased_coords, [], 2));
matlab_x = round(rebased_coords(1,:) + trans(1) + 1);
matlab_y = round(rebased_coords(2,:) + trans(2) + 1);


% Create 2D array for convolution of lattice points, compute 1D indices, fill array
L = zeros(max(matlab_x), max(matlab_y));
matlab_i = sub2ind(size(L), matlab_x, matlab_y);
L(matlab_i) = z;


% Compute output
C = conv2(prefilt, L);
coeffs = C(:);


% Compute coordinates for output lattice points in original basis
prefilt_trans = (size(prefilt) - 1) ./ 2;

xi = 1:size(C,1);
yi = 1:size(C,2);
[XI,YI] = ndgrid(xi - trans(1) - prefilt_trans(1) - 1, yi - trans(2) - prefilt_trans(2) - 1);

output_coords = bases * [XI(:)'; YI(:)'];
cx = output_coords(1,:);
cy = output_coords(2,:);