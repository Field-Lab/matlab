function kx = filter_stimulus_direct(K,X,basepars,dt)

kx = zeros(size(X,2),basepars.Nneurons);
1;
for n=1:basepars.Nneurons
    
    if (iscell(basepars.crop_idx))
        crop_idx_n = basepars.crop_idx{n};
    else
        crop_idx_n = basepars.crop_idx(:,n);
    end
    %kx(:,n) = sum(fastconv(X(crop_idx_n,:),K(:,:,n),basepars.n,size(X,2)),1);
1;
    kx(:,n) = dt .* sum(fastconv(X(crop_idx_n,:),K(:,:,n),length(crop_idx_n),size(X,2),basepars.padval),1);
end

1;