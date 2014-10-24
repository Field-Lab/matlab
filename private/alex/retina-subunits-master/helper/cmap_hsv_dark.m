function c = cmap_hsv_dark(N)

H = linspace(0,1,N+1);
H = H(1:end-1);
%H = linspace(0,1,N+2);
S = ones(1,N)*0.6;
V = ones(1,N)*0.8;
c = hsv2rgb([H(:) S(:) V(:)]);
%c = c(1:end-2,:);
