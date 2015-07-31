function processGuiOutput(currVec,iter,ratios,time, randflg, gapParam, saveFile, saveLocation)
	%PL = Pulse Library File
	%PLI = Pulse Library Index File
	%ESF = Event Sequence File
	%Make sure currVec is in micro-amperes (uA)

	function out = daqify(vec,range,numBits) %range is fixed for whole PL for now
		%http://hyperphysics.phy-astr.gsu.edu/hbase/electronic/dac.html
		numVals = 2^numBits;
		lvec = length(vec);
		ampVals = linspace(0,range,numVals); %assuming littlest bit corresponds to highest value
		out = zeros(lvec,numBits+1); %extra bit for polarity
		for i=1:lvec
			s = sign(vec(i));
			if s == -1; s = 1; else; s = 0; end %get polarity bit (not sure if 1 means negative)
			idx = closest(abs(vec(i)), ampVals);
			grot = fliplr(de2bi(idx-1));
			out(i,:) = [s padarray(grot,[0 (numBits-length(grot))],'pre')];
		end

		function c = closest(val,vec)
			[junk idx] = min(abs(vec - val));
			c = idx; 
		end
	end

	function plVecInt = make_plVecInt(ampVec, nratios, ntime)
		plVecInt = [];
		for i=ampVec
			for j=1:numel(nratios)
				for k=1:ntime(j)
					plVecInt = [plVecInt nratios(j)*i];
				end
			end
		end
	end

	function [ranges, pl_addition] = make_ranges(ampVec,rangeVec)
		grot = rangeVec(:,2);
		ranges = zeros(length(ampVec),2);
		for i=1:length(ampVec)
			val = find(diff(sign(grot - ampVec(i)))); %smallest value greater than 
			ranges(i,:) = rangeVec(val(1),:); %in case there is an equivalence and you get two values from the line above
		end
		pl_addition = 0; 
%		uranges = unique(ranges,'rows');
%		for i=1:length(uranges)
	end	


	%Parameters
	timeRes = 50; %microseconds
	numChannel = 512;
	numBits = 7;
	preDaqVec = [0 0 0 0 1 0 0 0]; %bits other than daq values for stimulation (only 'connect')
	range = 4; %microamps. Range is ASSUMED to be fixed for now
	extraCommandNum = 0; %commands other than stimulation pulses
	rangeVec = [-1 .06; -2 .250; -3 1; -4 4; -5 16; -6 64; -7 250; -8 1000]; %right column is in microamps

	%Make amplitude Vector
	ampVec = linspace(currVec(1),currVec(3),floor((currVec(3) - currVec(1))/currVec(2))+1); 

	%Normalize
	nratios = ratios/max(abs(ratios)); 		
	ntime = time/timeRes; 
	pulseLen = sum(time); %assume fixed pulse length for all pulses

	%Make whole number version of PL
	plVecInt = make_plVecInt(ampVec, nratios, ntime);
	[ranges,pl_additon] = make_ranges(ampVec,rangeVec);
	%Make PL
	plVec = daqify(plVecInt, range,numBits);
	plVec = [repmat(preDaqVec, size(plVec,1), 1) plVec];
	pl = typecast(uint16(bi2de(plVec)),'int16');

	%Make PLI
	pli = [1:extraCommandNum linspace(0,numel(plVecInt),numel(ampVec)+1)+1+extraCommandNum]; %assume each command is only 50 us
	pli = int32(pli(1:numel(pli)-1)); %discard last element

	%Make ESF
	chanVec = 1:numChannel;
	esfmat = zeros(iter*numel(ampVec)*numChannel + extraCommandNum*numChannel, 3);
	%random scan
	if randflg == 1
		esfmatInd = 1;
		timInc = 0;
		for i=1:length(pli) %same as numel(ampVec)
			for j=1:iter
				nchanVec = chanVec(randperm(length(chanVec)));
				for k=1:numChannel
					esfmat(esfmatInd,:) = [timInc nchanVec(k) i];
					esfmatInd = esfmatInd + 1;
					timInc = timInc + pulseLen + gapParam;
				end
			end
		end
	end
	
	%save everything
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
