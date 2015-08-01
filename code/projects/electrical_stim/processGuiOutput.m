function processGuiOutput(currVec,iter,ratios,time, randflg, gapParam, saveFile, saveLocation)
	%PL = Pulse Library File
	%PLI = Pulse Library Index File
	%ESF = Event Sequence File
	%Make sure currVec is in micro-amperes (uA)
	%Parts in this code that make assumptions that output is random noise:
   		%	

	function out = daqify(vec, ranges, numBits, plVecInd) 
		%http://hyperphysics.phy-astr.gsu.edu/hbase/electronic/dac.html
		numVals = 2^numBits;
		lvec = length(vec);
		out = zeros(lvec,numBits+1); %extra bit for polarity
		for i=1:lvec
			ampVals = linspace(0,ranges(plVecInd(i)),numVals); %assuming littlest bit corresponds to highest value
			s = sign(vec(i));
			if s == -1; s = 1; else s = 0; end %get polarity bit (not sure if 1 means negative)
			idx = closest(abs(vec(i)), ampVals);
			grot = fliplr(de2bi(idx-1));
			out(i,:) = [s padarray(grot,[0 (numBits-length(grot))],'pre')];
		end

		function c = closest(val,vec)
			[~,idx] = min(abs(vec - val));
			c = idx; 
		end
	end

	function [plVecInt,plVecInd] = make_plVecInt(ampVec, nratios, ntime)
		plVecInt = [];
		plVecInd = [];
		for i=1:length(ampVec)
			for j=1:numel(nratios)
				for k=1:ntime(j)
					plVecInt = [plVecInt nratios(j)*ampVec(i)];
					plVecInd = [plVecInd i];
				end
			end
		end
	end

	function ranges = make_ranges(ampVec,rangeVec)
		grot = rangeVec(:,2);
		ranges = zeros(length(ampVec),2);
		for i=1:length(ampVec)
			val = find(diff(sign(grot - ampVec(i)))); %smallest value greater than 
			ranges(i,:) = rangeVec(val(1)+1,:); %in case there is an equivalence and you get two values from the line above
		end
	end	

	function [esfmat,finalTimInc] = make_esfmat_random(numChannel, ampVec, iter, extraCommandNum, ranges, gapParam, pulseLen, simulCmd)
		chanVec = 1:numChannel;
		esfmat = zeros(iter*numel(ampVec)*numChannel + extraCommandNum*numChannel, 3);
		esfmatInd = 1;
		timInc = 0;
		for y=1:length(ampVec) 
			%set range in all channels
			rangeNum = ranges(y,1);
			[esfmat, timInc, esfmatInd] = set_range_all_channels(rangeNum, esfmat, esfmatInd, timInc, numChannel, simulCmd);
			for x=1:iter
				nchanVec = chanVec(randperm(length(chanVec)));
				for k=1:numChannel
					esfmat(esfmatInd,:) = [timInc nchanVec(k) y];
					esfmatInd = esfmatInd + 1;
					timInc = timInc + pulseLen + gapParam;
				end
			end
		end
		finalTimInc = timInc;

		function [esfmat, timInc, esfmatInd] = set_range_all_channels(rangeNum, esfmat, esfmatInd, timInc, numChannel, simulCmd)
			chanVec = 1:numChannel;
			chanInd = 1;
			for i=1:(numChannel/simulCmd)
				for j=1:simulCmd		
					esfmat(esfmatInd,:) = [timInc chanVec(chanInd) rangeNum];
					chanInd = chanInd + 1;
					esfmatInd = esfmatInd + 1;
				end
				timInc = timInc + timeRes;
			end
		end
	end

	function adjust_esfmat_times(esfmat, timeBreak, finalTimInc)
		esfTimes = esfmat(:,1);
		lEsfT = length(esfTimes);
		maxTime = esfTimes(end);
		lastTimeMult = ceil(maxTime/timeBreak);
		breakVec = linspace(timeBreak, timeBreak*lastTimeMult, lastTimeMult);
		for i=breakVec(1:(end-1))
			tind = find((esfTimes - i) < 0, 1, 'last');
			if esfTimes(tind) ~= i %if time doesn't magically correspond to a break, then shift everything downwards
				tgap = i - esfTimes(tind);
				for identical_inds = find(esfTimes==esfTimes(tind)) %first deal with all time stamps identical to the current one
					esfTimes(identical_inds) = esfTimes(identical_inds)+tgap;
				end
				for inds = (tind+1):lEsfT %deal with all inds after this one
					esfTimes(inds) = esfTimes(inds)+tgap;
				end
				if esfTimes(end) > timeBreak*lastTimeMult
			end
		end
		%deal with last value separately
		if finalTimInc



	%Parameters--------------------------------------------------
	timeRes = 50; %microseconds
	timeBreak = 500000; %microseconds that should not be crossed (no event should be occuring through increments of this time)
	simulCmd = 8;
	numChannel = 512;
	numBits = 7;
	preDaqVec = [0 0 0 0 1 0 0 0]; %bits other than daq values for stimulation (only 'connect')
	%range = 4; %microamps. Range is ASSUMED to be fixed for now
	extraCommandNum = 0; %commands other than stimulation pulses
	rangeVec = [-1 .06; -2 .250; -3 1; -4 4; -5 16; -6 64; -7 250; -8 1000]; %right column is in microamps
	%Make amplitude Vector
	ampVec = linspace(currVec(1),currVec(3),floor((currVec(3) - currVec(1))/currVec(2))+1); 
	%Normalize
	nratios = ratios/max(abs(ratios)); 		
	ntime = time/timeRes; 
	pulseLen = sum(time); %assume fixed pulse length for all pulses
	%paramCell = {timeRes numChannel numBits extraCommandNum pulseLen ampVec gapParam}; %for use in ESF

	%Make PL--------------------------------------------------
	%Make whole number version of PL
	[plVecInt,plVecInd] = make_plVecInt(ampVec, nratios, ntime);
	ranges = make_ranges(ampVec,rangeVec);
	%Make PL
	plVec = daqify(plVecInt, ranges, numBits, plVecInd);
	plVec = [repmat(preDaqVec, size(plVec,1), 1) plVec];
	pl = typecast(uint16(bi2de(plVec)),'int16');

	%Make PLI--------------------------------------------------
	pli = [1:extraCommandNum linspace(0,numel(plVecInt),numel(ampVec)+1)+1+extraCommandNum]; %assume each command is only 50 us
	pli = int32(pli(1:numel(pli)-1)); %discard last element

	%Make ESF--------------------------------------------------
	if randflg
        esfmat = make_esfmat_random(numChannel, ampVec, iter, extraCommandNum, ranges, gapParam, pulseLen, simulCmd); 
		esfmat = adjust_esfmat_times(esfmat, timeBreak, finalTimInc);
        disp(esfmat)
    end
	
	%Save everything--------------------------------------------------
	fid = fopen([saveLocation filesep saveFile '.slf'], 'w', 'l');
	fwrite(fid, pl, 'int16');
	fclose(fid);

	fid = fopen([saveLocation filesep saveFile '.sif'], 'w', 'l');
	fwrite(fid, pli, 'int32');
	fclose(fid);

	fid = fopen([saveLocation filesep saveFile '.sef'], 'w', 'l');
	fwrite(fid, esfmat, 'int32');
	fclose(fid);
end
