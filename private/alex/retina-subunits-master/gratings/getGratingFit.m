function [gfit fit] = getGratingFit(g,fit,opts,refit)

if ~exist('refit','var');
	refit = 0;
end

% find best fitting scale and shift 

shift = [0:29];
mse = zeros(length(shift),1);
for ishift=1:length(shift)
	gTmp = shiftGratings(g,shift(ishift));
	dataTmp = getGratingData(gTmp,opts);
	out = evalFit(dataTmp,fit);
	mse(ishift) = sum((out.Z_t*0.15-dataTmp.R_t).^2);
end
[~, mnInd] = min(mse);
shift_hat = shift(mnInd);

% apply the best shifting
gTmp = shiftGratings(g,shift_hat);

data = getGratingData(gTmp,opts);
out = evalFit(data,fit);
b = fitLine(out.Z_t,data.R_t);

% reformat data and refit nonlinearities
if refit
	data = getGratingData(gTmp,opts);
	switch fit.model
		case 'subunit'
			fit.error = 'mse';
			for i=1:1
				fit = fitF(data,fit);
				%fit = fitG(data,fit);
			end
		case 'ln'
			fit.error = 'mse';
			%fit = fitG(data,fit);
	end
end

% evalute fit on subset
opts.cycles = 1;
data = getGratingData(gTmp,opts);
reliability = repeatReliability(data.R_t);
out = evalFit(data,fit);
cr_boot_all = zeros(1,size(data.R_t,1));
for i=1:size(data.R_t,1)
	tmp = corrcoef(out.Z_t,data.R_t(i,:));
	cr_boot_all(i) = tmp(1,2)^2;
end
cr_boot = nanmean(cr_boot_all);
cr = cr_boot / reliability;


% evaluate fit on full data set (for plotting)
data = getGratingData(gTmp);
out = evalFit(data,fit);
if ~refit
	out.Z_t = out.Z_t*b(1) + b(2);
	gfit.b = b;
end

% collect results
gfit.stim = gTmp.stim;
gfit.shift_hat = shift_hat;
gfit.cycleAve = reshape(out.Z_t,[g.cycleDur g.nPer g.nPh]);
gfit.cr = cr;


