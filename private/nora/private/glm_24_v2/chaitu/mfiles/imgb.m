function imgb(A, norm, t, dims)
% function imgb(A, norm, transpose, dims)
%
%     A:     matrix of column vectorized images
%     norm:  stretch each to fill range (but keep 0 at 0) or not?
%     transpose:  transpose x and y of each image
%     dims:  [rows cols] of the resulting image arrangement
% 
if nargin < 2, norm = 1; end;

if nargin < 3, t = 0; end;

if nargin < 4, dims = []; end;

if ndims(A)==3 & size(A,3)>1,
  A = reshape(permute(A,[1 3 2]),[size(A,1),size(A,2)*size(A,3)]);
end;

imagesc(basis2img(A,norm,dims,t)); axis image off;
set(gca,'clim',[-1 1]);


function I = basis2img(A, normalized, dims, tr);
% function I = basis2img(A, normalized, dims, transpose);
% 
%    normalized   0: no rescaling of bfs
%                 1: stretch-to-fill each bfs but preserve 0
%                 [min max]: rescale global axes to given range
%
%    dims         dimensions row/col of output array
%
%    tr           flag to transpose each basis function
%
%
%    I            an image of 2D bfs

N = size(A,1);
while (sqrt(N)-floor(sqrt(N))),
  A = [A; zeros(1,size(A,2))];
  N = size(A,1);
end;
M = size(A,2);
w = sqrt(N);

if nargin < 2 | isempty(normalized),
  normalized = 1;
end;

if nargin < 4, tr = 0; end;

if length(normalized)>1,
  rangemin = normalized(1);
  rangemax = normalized(2);
  normalized = 0;
else,
  rangemin = min([min(A(:)) -max(A(:))]); 
  rangemax = -rangemin;
end;
  
spacing = 1;
if nargin < 3 | isempty(dims),
  a = floor(sqrt(size(A,2)));
  if size(A,2) > a*(a+1),  nrows = a+1;
  else,                    nrows = a; end;
  ncols = ceil(sqrt(size(A,2)));
else,
  nrows = dims(1); ncols = dims(2); 
end;
  

sI = [w*nrows+(nrows-1)*spacing w*ncols+(ncols-1)*spacing];
I = -1*ones(sI(1), sI(2));



if ~normalized, % stretch all together  
  A = (A-rangemin)./(rangemax-rangemin); 
  A(find(A<0))=0; A(find(A>1))=1;
  A = 2*A-1;
end;

for i = 1:size(A,2),
  patch = reshape(A(:,i), w, w);
  if normalized, % stretch each patch
    prange = max(abs(patch(:)));
    if prange > .00001, patch = patch / prange; end;
  end;

  sx = (w+spacing) * rem(i-1,ncols)+1;
  sy = (w+spacing) * floor((i-1)/ncols)+1;
  if tr, patch = patch'; end;
  I(sy:sy+w-1, sx:sx+w-1,:) = patch;
end

if nrows*ncols > M,
  ix1 = (nrows-1)*(w+1)+1:sI(1);
  ix2 = mod(M,ncols)*(w+1)+1:sI(2);
  msk = -.5*ones(length(ix1), length(ix2)); 
  msk(1:2:end,1:2:end) = -1;
  msk(2:2:end,2:2:end) = -1;
  I(ix1,ix2) = msk;
end;
