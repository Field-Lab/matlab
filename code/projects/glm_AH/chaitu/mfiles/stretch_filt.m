function K2 = stretch_filt(K,r)

    [n m] = size(K);
    K2 = zeros(r^2*n,m);
    
    for j=1:size(K2,2)
        
        K2(:,j) = reshape(imresize(reshape(K(:,j),sqrt(n),sqrt(n)),r,'nearest'), r^2 * n,1);
        
    end
    
    
    K2 = K2 .* norm(K,'fro')/norm(K2,'fro');