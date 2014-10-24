function getTimeScaling(celldat,train,test,fit)


	[simdat simtrain simtest] = getSimulation(celldat,train,test,fit,size(train.X_ct,2),size(test.X_ct,2));

	
	coneInputs = zscore(dat.coneInputs);
	coneInputs = tempFilterCones(coneInputs,celldat.timeCourse/norm(celldat.timeCourse));

	[fulldattrain.X_ct fulldattest.X_ct] = partitionDat(coneInputs,celldat.trainFrac,celldat.subDivide,1);

	fulldattrain.X_ct = fulldattrain.X_ct';
	fulldattest.X_ct = fulldattest.X_ct';


	count = 1;
	for i=round(linspace(10000,60000,10));
		for rep=1:5
			rng = [1:i];
			tmp = bsxfun(@times,fulldattrain.R_t(:,simtrain.R_t(rng)>0),simtrain.Z_t(rng(simtrain.R_t(rng)>0)));
			tmp = fulldattrain.X_ct(:,simtrain.R_t(rng)>0);
			sta = mean(tmp,2);

			rng = [1:i];
			tmp2 = bsxfun(@times,fulldattrain.X_ct(:,train.R_t(rng)>0),train.R_t(rng(train.R_t(rng)>0)));
			tmp2 = fulldattrain.X_ct(:,train.R_t(rng)>0);
			sta2 = mean(tmp2,2);

			snr1(count,rep) = max(abs(sta))/(1.4826*mad(sta,1))
			snr2(count,rep) = max(abs(sta2))/(1.4826*mad(sta2,1))
		end
		count = count + 1;
	end