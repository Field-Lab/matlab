function vects_norm = normalize_to_unit_sphere(vects)
% normalize a vertical list of vectors to each  lie on the unit sphere
%
% gauthier 2008-10
%

% compute normalization factors
norm_factors = sqrt(sum(vects.^2,2));

% divide 
vects_norm = vects./repmat(norm_factors,1,size(vects,2));
