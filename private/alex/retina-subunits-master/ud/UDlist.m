

runids = {'2011-10-14-1-data001-0','2011-10-25-5-data001-0','2012-04-13-1-data002',...
'2011-10-25-9-data006-0','2012-04-13-1-data002','2012-09-06-0-data004','2012-09-13-2-data005'};
rgcids = {[136,947,6811,7471],[184,454,698,948,2554,2583,3169,6408],[3002],[4066,4489,5436],...
[2389, 6286, 3002],[7008],[6991]};

respAll = [];
predicAll = [];
rsing = [];
subflag = [];
for i=1:length(runids)
	tmpids = rgcids{i};
	for j=1:length(tmpids)
		[histx histy resSUB resSING coneids] = getUDfiles(runids{i},tmpids(j));
		if resSUB.fit_SUB.S ~= resSUB.fit_SUB.C
			subflag = [subflag, 1];
		else
			subflag = [subflag, 0];
		end
		[resp predic] = getUD(histx,histy,coneids,resSING.fit_SUB,4);
		resp = resp(~logical(diag(ones(1,4))));
		predic = predic(~logical(diag(ones(1,4))));
		respAll = [respAll, resp(:)];
		predicAll = [predicAll, predic(:)];
		foo = corrcoef(resp,predic);
		rsing = [rsing;foo(1,2)]
	end
end

predicAllVecSub = predicAll(:);
respAllVecSub = respAll(:);

% for n = 21, R is 0.48 for subunit, 0.32 for LN, and 0.32 for Singletons
predicAllVecLN = predicAll(:);
respAllVecLN = respAll(:);
n = length(predicAllVecSub);
for iboot=1:1000
	bootsmp = ceil(rand(1,n)*n);
	foo1 = predicAllVecLN(bootsmp);
	foo2 = respAllVecLN(bootsmp);
	[r1 p] = corrcoef(foo1,foo2);
	foo3 = predicAllVecSub(bootsmp);
	foo4 = respAllVecSub(bootsmp);
	[r2 p] = corrcoef(foo3,foo4);
	diffboot(iboot) = r2(1,2) - r1(1,2);
end
