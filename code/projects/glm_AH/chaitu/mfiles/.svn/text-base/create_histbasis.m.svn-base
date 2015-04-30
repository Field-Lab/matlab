function [basis,t_psb] = create_histbasis(alpha,beta,bstretch,nofilters,Mhist,dt,norm_mode,spacing)

%basis =  construct_cosine_basis2(alpha,beta,bstretch,nofilters,dt,spacing); % raised cosine bumps
[basis,~,t_psb] =  construct_cosine_basis3(alpha,beta,bstretch,nofilters,dt,spacing); % edoi, 2011-01-07.

if (alpha < 0)
    offset = floor(abs(alpha)/dt);
    basis1 = basis(offset+1:min(size(basis,1),offset+Mhist),:);
    basis = [basis1; zeros(max(0,Mhist-size(basis1,1)),size(basis,2))];
else
    basis = [zeros((Mhist-size(basis,1)),size(basis,2)); basis];
end


switch(norm_mode)
    case 'area'
        areas = sum(basis,1); % Normalization by signed area
    case 'L2'
        areas = sqrt(sum(basis.^2,1)); % Normalization by L2 norm
    case 'none'
        areas = ones(1,size(basis,2));
end
basis = basis ./ repmat(areas,Mhist,1); % normalize the basis

basis(isnan(basis)) = 0;