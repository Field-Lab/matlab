function data = getGratingData(g,opts)

	if ~exist('opts','var')
		opts = struct;
	end

	if ~isfield(opts,'cycles')
		opts.cycles = 0;
	end

	if ~isfield(opts,'pers')
		opts.pers = [1:g.nPer];
	end

	%foo = repmat(g.stim,[g.nCycle 1 1 1]); if we want to work with g.timeCourses...

	% do we want single cycles or cycle averages
	if opts.cycles
			data.R_t = zeros(g.nCycle,g.cycleDur*length(opts.pers)*g.nPh);
		for icycle = 1:g.nCycle
			rng = (icycle-1)*g.cycleDur+1:(icycle)*g.cycleDur;
			tmp = squeeze(mean(g.timeCourses(rng,:,opts.pers,:),2));
			data.R_t(icycle,:) = reshape(tmp,[g.cycleDur*length(opts.pers)*g.nPh 1])';
		end
	else
		data.R_t = reshape(g.cycleAve(:,opts.pers,:,:),[g.cycleDur*length(opts.pers)*g.nPh 1])';
	end

data.X_ct = reshape(g.stim(:,opts.pers,:,:),[g.cycleDur*length(opts.pers)*g.nPh size(g.stim,4)])';