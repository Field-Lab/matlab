% Return an orthogonal basis that has the same span as the orogiinal. (uses
% QR factorization)

function obasis = orthogonalize_basis(basis)


[Q R] = qr(basis,0);
obasis = Q(:,1:size(basis,2));