function [coeffs2 obasis] = get_ortho_coeffs(coeffs1,basis1)

obasis = orthogonalize_basis(basis1);
coeffs2 = (obasis'*obasis)\(obasis'* (basis1*coeffs1));
1;

if (norm(obasis*coeffs2 - basis1*coeffs1)>1e-12)
    error('orthogonal coeffs are inaccurate. discrepancy is %0.3f',norm(obasis*coeffs2 - coeffs1*basis1));
end