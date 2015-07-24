function [csta cidx] = get_cropped_sta(sta,height,width,sidelen)


zmeanSTA = reshape(sum((sta - repmat(mean(sta,2),1,size(sta,2))).^2,2),height,width);
1;
cidx = get_crop_idx(zmeanSTA,sidelen); % get the relevant indices
csta = sta(cidx,:);
