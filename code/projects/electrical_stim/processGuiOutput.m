function processGuiOutput(currVec,iter,ratios,time, randflg, randampflg, gapParam, saveFile, saveLocation)
	%PL = Pulse Library File
	%PLI = Pulse Library Index File
	%ESF = Event Sequence File
	%Make sure currVec is in micro-amperes (uA)
	%Parts in this code that make assumptions that output is random noise:
   		%	

	function out = daqify(vec, ranges, numBits, plVecInd, ist_len) 
		%http://hyperphysics.phy-astr.gsu.edu/hbase/electronic/dac.html
		plVecInd = [plVecInd 1e9]; %for indexing purposes
		numVals = 2^numBits;
		lvec = length(vec);
		out = zeros(lvec,ist_len); %extra bit for polarity
		con_pre = [0 0 0 0 1 0 0 0]; 
		none_pre = [0 0 0 0 0 0 0 0]; 
		prev_ind = 0;
		for i=1:lvec
			ampVals = linspace(0,ranges(plVecInd(i),2),numVals); %assuming littlest bit corresponds to highest value
			s = sign(vec(i));
			if s == -1; s = 0; else s = 1; end %get polarity bit 
			[idx, val] = closest(abs(vec(i)), ampVals);
			%disp([num2str(vec(i)) ' gets ' num2str(val) ' and ' num2str(idx) ' for ' num2str(ampVals)])
			grot = fliplr(de2bi(idx-1));
			%decide whether to connect or not
			if (plVecInd(i) ~= prev_ind) || (plVecInd(i) ~= plVecInd(i+1))
				out(i,:) = [none_pre s padarray(grot,[0 (numBits-length(grot))],'pre')];
			else
				out(i,:) = [con_pre s padarray(grot,[0 (numBits-length(grot))],'pre')];
			end
			prev_ind = plVecInd(i);
		end

		function [idx,val] = closest(val,vec)
			[val,idx] = min(abs(vec - val));
		end
	end

	function [plVecInt,plVecInd] = make_plVecInt(ampVec, nratios, ntime)
		plVecInt = [];
		plVecInd = [];
		for i=1:length(ampVec)
			for j=1:numel(nratios) %to account for flanking regions
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

	function [esfmat,finalTimInc, esfmat_type] = make_esfmat_random(numChannel, ampVec, iter, ranges, gapParam, pulseLen, frames_for_all_channels, randampflg)
		extraCommandNum = numel(unique(ranges(:,2)));
		chanVec = 1:numChannel;
		esfmat = zeros(iter*numel(ampVec)*numChannel + extraCommandNum, 3);
		esfmat_type = zeros(iter*numel(ampVec)*numChannel + extraCommandNum, 1);
		esfmatInd = 1;
		timInc = 1;

		prev_rangeNum = 0;
		%Eventually have to make this have random amplitude, in which case
		%we would change the range of the channels individually, not all at once
		%this is inefficient by least a factor of 8 times the proportion of time
		%spent on changing ranges, so will implement this later.
		for y=1:length(ampVec) 
			%set range in all channels
			rangeNum = ranges(y,1);
			if rangeNum ~= prev_rangeNum
				esfmat(esfmatInd,:) = [timInc 0 rangeNum];
				esfmat_type(esfmatInd) = 1;
				esfmatInd = esfmatInd + 1;
				timInc = timInc + frames_for_all_channels;
			end
			prev_rangeNum = rangeNum;

			for x=1:iter
				nchanVec = chanVec(randperm(length(chanVec)));
				for k=1:numChannel
					esfmat(esfmatInd,:) = [timInc nchanVec(k) y];
					esfmat_type(esfmatInd) = 2;
					esfmatInd = esfmatInd + 1;
					timInc = timInc + pulseLen + gapParam;
				end
			end
		end
		finalTimInc = timInc;
	end

	function esfmat = adjust_esfmat_times(esfmat, esfmat_type, timeBreak, pulseLen, frames_for_all_channels)
		lEsfT = length(esfmat);
		esfTimes = [esfmat(:,1) esfmat(:,1)];
		for i = 1:lEsfT
			if esfmat_type(i) == 1; esfTimes(i,2) = esfTimes(i,2) + frames_for_all_channels;
			else esfTimes(i,2) = esfTimes(i,2) + pulseLen; end;
		end
		breakVal = timeBreak;
		while esfTimes(end) > breakVal
			for i = 1:length(esfTimes)
				if (esfTimes(i,1) < breakVal) && (esfTimes(i,2) > breakVal)
					tgap = breakVal - esfTimes(i,1);
					esfTimes(i:end,:) = esfTimes(i:end,:) + tgap;
				end
			end
			breakVal = breakVal + timeBreak; 
		end
	end
		




	%Parameters--------------------------------------------------
	ist_len = 16;
	timeRes = 50; %microseconds. Also taken to be command time
	gapParam = gapParam/timeRes; %convert to frames units
	%timeBreak = 500000; %microseconds that should not be crossed (no event should be occuring through increments of this time)
	timeBreak = 10000; %frames that should not be crossed (no event should be occuring through increments of this time)
	simulCmd = 8;
	numChannel = 512;
	frames_for_all_channels = numChannel/simulCmd;
	numBits = 7;
	preDaqVec = [0 0 0 0 1 0 0 0]; %bits other than daq values for stimulation (only 'connect')
	%range = 4; %microamps. Range is ASSUMED to be fixed for now
	rangeVec = [-1 .06; -2 .250; -3 1; -4 4; -5 16; -6 64; -7 250; -8 1000]; %right column is in microamps
	%Make amplitude Vector
	ampVec = linspace(currVec(1),currVec(3),floor((currVec(3) - currVec(1))/currVec(2))+1); 
	%Normalize
	nratios = ratios/max(abs(ratios)); nratios = [nratios(1) nratios nratios(end)]; %to account for flanking regions		
	ntime = time/timeRes; ntime = [1 ntime 1]; %to account for flanking regions
	%pulseLen = sum(time); %assume fixed pulse length for all pulses
	pulseLen = sum(ntime); %assume fixed pulse length for all pulses

	%Make PL--------------------------------------------------
	%Make whole number version of PL
	[plVecInt,plVecInd] = make_plVecInt(ampVec, nratios, ntime);
	ranges = make_ranges(ampVec,rangeVec);
	%Make PL
	plVec = daqify(plVecInt, ranges, numBits, plVecInd, ist_len);
	plVec = [repmat(preDaqVec, size(plVec,1), 1) plVec];
	pl = typecast(uint16(bi2de(plVec)),'int16');

	%Make PLI--------------------------------------------------
	pli = [linspace(0,numel(plVecInt),numel(ampVec)+1)+1]; %assume each command is only 50 us
	pli = int32(pli(1:numel(pli)-1)); %discard last element

	%Make ESF--------------------------------------------------
	if randflg
        [esfmat, finalTimInc, esfmat_type] = make_esfmat_random(numChannel, ampVec, iter, ranges, gapParam, pulseLen,frames_for_all_channels, randampflg); 
		esfmat = adjust_esfmat_times(esfmat, esfmat_type,timeBreak, pulseLen, frames_for_all_channels);
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
