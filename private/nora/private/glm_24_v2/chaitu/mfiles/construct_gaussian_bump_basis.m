% Function to construct a set of 2D gaussian basis functions  with std. devs equally spaced in [min_std,max_std]


% Arguments:

% centers - will create a set of basis functions centered around each point in centers 
%         - should be a n x 2 matrix with row/col indices (in the s x s box) in the 1st/2nd col, resp. 
% min_std,max_std: minimum and maximum standard deviation
% nfilters - number of filters PER center
% s - length of one side of the square cropped box
% plotFlag - whether or not to plot the basis

function F = construct_gaussian_bump_basis(centers,min_std,max_std,nfilters,s,plotFlag)

if (isempty(centers))
    centers = [floor(s/2)+1 floor(s/2)+1];
end

stds = linspace(min_std,max_std,nfilters);
F = zeros(s^2,nfilters*size(centers,1));



eff_len = 2*s+1;

%tl_row = s/2;
%tl_col = s/2;

c_row = s+1;
c_col = s+1;

counter = 1;
for i=1:size(centers,1)
    
    for j=1:nfilters
        
         Fjraw = fspecial('gaussian',eff_len,stds(j));

        1;
         % Crop the window based on center
         Fj = Fjraw(c_row-centers(i,1)+1:c_row-centers(i,1)+s,c_col-centers(i,2)+1:c_col-centers(i,2)+s);
         %Fj = Fjraw(tl_row+centers(i,1)-floor(s/2):tl_row+centers(i,1)+floor(s/2)-1,tl_col+centers(i,2)-floor(s/2):tl_col+centers(i,2)+floor(s/2)-1);
         F(:,counter) = reshape(Fj,s^2,1);
        
        if (nargin > 5 && plotFlag)
            subplot(size(centers,1),nfilters,(i-1)*nfilters+j), imagesc(Fj), axis image;
        end
        
        F(:,counter) = F(:,counter)./norm(F(:,counter)); % normalize by L2 norm
        counter = counter + 1;
        
    end
end
F(isnan(F)) = 0;