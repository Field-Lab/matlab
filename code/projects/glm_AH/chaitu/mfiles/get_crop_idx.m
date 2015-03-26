function [idx box_m box_n] = get_crop_idx(A,sidelen,forceFlag)

% Given a 2D nonneg-valued matrix, we want return the indices of a box around the maximum point.
% The box has dimensions sidelen x sidelen.
% If the forceFlag is on, then if the maximum point is on the border of the image, we return a truncated cropped index.


if (~exist('forceFlag','var'))
    forceFlag = 0;
end

%[h w] = size(A); % note: stimulus in snle is conventionally defined w x h (edoi)
[I, J] = size(A); 

%Avec = A(:);
%[~, sorted_idx] = sort(Avec,'descend');
[~, ma_idx] = max((abs(A(:))));

[mi, mj] = ind2sub([I,J], ma_idx); % indices of the max elem

lb = sidelen/2;
ub = (sidelen/2 - 1);

if (mi > lb && mi < (I-ub) && mj > lb && mj < (J-ub))
    % The maximum is an interior point        
    rl = mi-floor(sidelen/2);
    rr = rl + sidelen-1;
    cl =  mj-floor(sidelen/2);
    cr = cl + sidelen-1;    
    rvec = rl:rr;
    cvec = cl:cr;
    fprintf('Maximum is at index (%d,%d). Returning the box between rows %d - %d and cols %d - %d\n',mi,mj,rl,rr,cl,cr);
    box_m = length(rvec);
    box_n = box_m;
    
    idx = sub2ind(size(A),repmat(rvec',length(cvec),1),reprows(cvec(:),length(rvec)));
    return;
end

if (forceFlag)
    % The maximum point is on the border of the image
    
    if (mi <= lb)
        rvec = 1:sidelen;
    elseif (mi > (I-ub))
        rvec = I-sidelen+1:I;
    else
        rvec = mi-floor(sidelen/2):mi-floor(sidelen/2)+sidelen-1;
    end
        
    if (mj <= lb)
        cvec = 1:sidelen;
    elseif (mj > (J-ub))
        cvec = J-sidelen+1:J;
    else
        cvec = mj-floor(sidelen/2):mj+floor(sidelen/2)+sidelen-1;
    end
    rveclen = length(rvec);
    cveclen = length(cvec);
    cvec = repcols(cvec,rveclen)';
    
    fprintf('Maximum is at index (%d,%d). Returning the TRUNCATED box between rows %d - %d and cols %d - %d\n',mi,mj,rvec(1),rvec(end),cvec(1),cvec(end));    
    1;
    idx = sub2ind(size(A),repmat(rvec',cveclen,1),cvec);
    
    box_m = rveclen;
    box_n = cveclen;
    
    return;
end

% Return an error
fprintf('ERROR: the receptive field of this neuron goes outside the stimulus boundary: mi=%d, mj=%d. stim size is %d x %d \n',mi,mj,I,J);
idx = [];

