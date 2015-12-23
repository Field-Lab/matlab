function esfmat_out = adjust_esfmat_times(esfmat, timeBreak, finalTimInc)
	esfTimes = [esfmat(:,1); finalTimInc];
	lEsfT = length(esfTimes);
	maxTime = esfTimes(end);
	lastTimeMult = ceil(maxTime/timeBreak);
	breakVal = timeBreak; %first time break
	while esfTimes(end) > breakVal
		tind = find((esfTimes - breakVal) <= 0, 1, 'last');
		if esfTimes(tind) ~= i %if time doesn't magically correspond to a break, then shift everything downwards
			tgap = breakVal - esfTimes(tind);
			%first deal with all time stamps identical to the current one
			identical_inds = find(esfTimes==esfTimes(tind)); 
			esfTimes(identical_inds) = esfTimes(identical_inds)+tgap;
			 %deal with all inds after this one
			remaining_inds = (tind+1):lEsfT;
			esfTimes(remaining_inds) = esfTimes(remaining_inds)+tgap;
		end
		breakVal = breakVal + timeBreak; 
		disp(esfTimes)
	end
	disp(esfTimes)
	esfmat(:,1) = esfTimes(1:end-1);
	esfmat_out = esfmat;
