function [space time svs] = get_spacetime_approx(A,r,plotflag,xdims)

if (size(A,1) == 1)
    space(1) = 1;
    time(:,1) = A;
    svs = 1;
    for i=2:r
        space(i) = 0;
        time(:,i) = 0;
        svs(i) = 0;
    end
    return;
end

if (size(A,2) == 1)
    space(:,1) = A;
    time(1) = 1;
    svs(1) = 1;
    for i=2:r
        space(i) = 0;
        time(:,i) = 0;
        svs(i) = 0;
    end    
    return;
end


[U S V] = svd(A);
n = size(A,1);
t = size(A,2);
1;
space = U(:,1:r);
time = V(:,1:r)*S(1:r,1:r);

for j=1:r
    [sorted sorted_idx] = sort(abs(space(:,j)),'descend');
    if (space(sorted_idx(1),j) < 0) % Flip signs to make space filter a positive bump
       space(:,j) = -space(:,j);
       time(:,j) = -time(:,j);
    end
end


svs = diag(S); svs = svs(1:r);
if (plotflag)
    
    if (~exist('xdims','var'))
        xdims = [sqrt(n) sqrt(n)];
    end
    %figure;
    for j=1:r
        subplot(2,r,j), imagesc(reshape(space(:,j),xdims)), axis image, colorbar;
        title(sprintf('Rank %f/%f',svs(j),sum(svs)));
        subplot(2,r,r+j), plot(time(:,j));
    end
end



