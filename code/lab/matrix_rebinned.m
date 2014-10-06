function rsta = matrix_rebinned(sta, bin_factor)
% matrix_rebinned     Bin neighboring stixels
%
%    This reduces the size of the matrix by bin_factor by summing the values of neighboring pixels
%
%
%
%
% usage:  sta = matrix_rebinned(sta, bin_factor)
%
% arguments:      sta - 4d matrix of an STA
%          bin_factor - integer
%
% outputs:        sta - rebinned STA
%
%




% BODY OF THE FUNCTION


% subsample the STA in each frame and color

for cc = 1:size(sta,3)
    for ff = 1:size(sta,4)
        
        

        %  this is done by summing over nxn regions

        M = sta(:,:,cc,ff);
        n = bin_factor;
        x_size = size(sta,2);
        y_size = size(sta,1);
        
        
        %figure(1);clf;imagesc(M)
        
        M = reshape(M,[n y_size/n x_size]);
        M = permute(M,[1 3 2]);
        M = reshape(M,[n*n x_size/n y_size/n]);
        M = squeeze(sum(M,1))';
        
        rsta(:,:,cc,ff) = M;
        
        %figure(2);clf;imagesc(M)
        
    end
end

