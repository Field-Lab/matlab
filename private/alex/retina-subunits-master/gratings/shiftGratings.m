function g = shiftGratings(g,shifter)

	for iper = 1:g.nPer
		for iph = 1:g.nPh
			g.stim(:,iper,iph,:) = circshift(squeeze(g.stim(:,iper,iph,:)),[-shifter 0]);
		end
	end

