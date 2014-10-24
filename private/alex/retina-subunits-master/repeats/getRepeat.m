function re = getRepeat(celldat,repeatrun)

% inputs are an rgcId, a dat structure for a WN run, a dat file for a repeatrun, 
% and a flag for whether we want to plot
% output is a struct with info about the repeat, and 

% find the cell
rgcId_re = repeatrun.cell_list_map{celldat.cellInd};

if isempty(rgcId_re)
	fprintf('(getRepeat) could not find matching cell for %g \n',celldat.rgcId);
	re = [];
	return;
end

cellInd_re = find(repeatrun.cell_ids==rgcId_re);

% get a matrix of repeat responses
resp = [];
nRepeats = length(repeatrun.spike_rate);
for i=1:nRepeats
	resp(i,:) = repeatrun.spike_rate{i}(cellInd_re,:);
end
reliability = repeatReliability(resp);
re.R_t = mean(resp);

% get cone inputs
re.X_ct = zscore(double(repeatrun.cone_inputs(:,celldat.coneIds)));
re.X_ct = tempFilterCones(re.X_ct,celldat.timeCourse/norm(celldat.timeCourse));
re.X_ct = re.X_ct';

% get model predictions
res = loadFit(celldat,'greedy',[]);

if isempty(res)
	fprintf('(getRepeat) could not find fit for %g \n',celldat.rgcId);
	re = [];
	return;
end

predic = evalFit(re,res.fit_SUB);
b = regress(re.R_t',[predic.Z_t; ones(1,length(predic.Z_t))]');
predic.Z_t = predic.Z_t*b(1) + b(2);
cr = corrcoef(predic.Z_t,re.R_t)^2;
cr = cr(1,2);

cr_boot_all = zeros(1,size(resp,1));
for i=1:size(resp,1)
	tmp = corrcoef(predic.Z_t,resp(i,:));
	cr_boot_all(i) = tmp(1,2)^2;
end
cr_boot = mean(cr_boot_all);

predic_LN = evalFit(re,res.fit_LN);
b = regress(re.R_t',[predic_LN.Z_t; ones(1,length(predic.Z_t))]');
predic_LN.Z_t = predic_LN.Z_t*b(1) + b(2);
cr_LN = corrcoef(predic_LN.Z_t,re.R_t)^2;
cr_LN = cr_LN(1,2);

cr_LN_boot_all = zeros(1,size(resp,1));
for i=1:size(resp,1)
	tmp = corrcoef(predic_LN.Z_t,resp(i,:));
	cr_LN_boot_all(i) = tmp(1,2)^2;
end
cr_LN_boot = mean(cr_LN_boot_all);


% collect outputs
re.reliability = reliability;
re.resp = resp;
re.predic = predic;
re.predic_LN = predic_LN;
re.cr = cr;
re.cr_LN = cr_LN;
re.cr_boot = cr_boot;
re.cr_LN_boot = cr_LN_boot;




