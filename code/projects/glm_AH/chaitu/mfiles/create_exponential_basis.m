function basis = create_exponential_basis(taus,Mhist,dt,norm_mode)


t = 0:dt:(Mhist-1)*dt;
basis = zeros(Mhist,length(taus));

for j=1:length(taus)
    basis(:,j) = exp(-t./taus(j));
end

switch(norm_mode)
    case 'area'
        areas = sum(basis,1); % Normalization by signed area
    case 'L2'
        areas = sqrt(sum(basis.^2,1)); % Normalization by L2 norm
end

basis = basis ./ repmat(areas,Mhist,1); % normalize the basis

basis(isnan(basis)) = 0;