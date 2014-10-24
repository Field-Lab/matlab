function out = repeatRel(resp)

% input is matrix of repeat responses, n repeats x time points
% output is the reliability estimates by correlating the average of
% all trials but 1 with the response on the remaining trial

for i=1:size(resp,1)
	inds_tmp = logical(zeros(1,size(resp,1)));
	inds_tmp(i) = 1;
	r_tmp = corrcoef(resp(inds_tmp,:),mean(resp(~inds_tmp,:)));
	r_all(i) = r_tmp(1,2)^2;
end

out = nanmean(r_all);
