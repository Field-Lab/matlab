function [resp predic] = getUD(histx,histy,coneids,fit,normOpt)

inds = histx > 0 & histx < 0.3;

for i=1:4
	for j=1:4
		resp(i,j) = mean(histy{i,j}(inds));
	end
end

for i=1:4
for j=1:4
	if i==j
		pairs{i,j} = coneids(i);
	else
		pairs{i,j}(1) = coneids(i);
		pairs{i,j}(2) = coneids(j);
	end
end
end

if isfield(fit,'k')
	nCones = length(fit.k);
else
	nCones = size(fit.locs_c,1);
end

stimAll = [];
respAll = [];
for i=1:4
for j=1:4
	stim = zeros(1,nCones);
	if length(pairs{i,j}) > 1 
		stim(pairs{i,j}(1)) = -4;
		stim(pairs{i,j}(2)) = 4;
		data.X_ct(:,1) = stim;
		data.R_t = 1;
		foo = evalFit(data,fit);
		predic(i,j) = foo.Z_t;
	else
		stim(pairs{i,j}) = -4;
		data.X_ct(:,1) = stim;
		data.R_t = 1;
		foo = evalFit(data,fit);
		predic(i,j) = foo.Z_t;
	end
	stimAll = [stimAll;stim]'';
	respAll = [respAll;resp(i,j)];
end
end

data.R_t = respAll';
data.X_ct = stimAll';

if normOpt == 1
	resp = resp / max(resp(:));
	predic = predic / max(predic(:));
end

if normOpt == 1
	resp = resp / max(resp(:));
	predic = predic / max(predic(:));
end
if normOpt == 2
	resp = zscore(resp,[],2);
	predic = zscore(predic,[],2);
end
if normOpt == 3
	resp = bsxfun(@minus,resp,min(resp,[],2));
	predic = bsxfun(@minus,predic,min(predic,[],2));
	resp = bsxfun(@rdivide,resp,max(resp,[],2));
	predic = bsxfun(@rdivide,predic,max(predic,[],2));
end
if normOpt == 4
	resp = bsxfun(@rdivide,resp,max(resp,[],2));
	predic = bsxfun(@rdivide,predic,max(predic,[],2));
end
if normOpt == 5
	resp = resp;
	fit = fitB(data,fit);
	foo = evalFit(data,fit);
	predic = reshape(foo.Z_t,4,4)';
end