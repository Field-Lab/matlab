function processGuiOutput(currVec,iter,ratios,time, randflg, gapParam, saveFile, saveLocation)
	%PL = Pulse Library File
	%PLI = Pulse Library Index File
	%ESF = Event Sequence File

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
			grot = fliplr(de2bi(idx));
			out(i,:) = [s padarray(grot,[0 (numBits-length(grot))],'pre')];
		end

		function c = closest(val,vec)
			[junk idx] = min(abs(vec - val));
			c = idx; 
		end
	end

	%Parameters
	timeRes = 50; %microseconds
	numChannel = 512;
	numBits = 7;
	setRecordVec = [0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0];
	preDaqVec = [0 0 0 0 1 0 0 0]; %bits other than daq values for stimulation (only 'connect')
	range = 1000; %microamps. Range is ASSUMED to be fixed for now
	extraCommandNum = 1; %commands other than stimulation pulses

	%Make amplitude Vector
	ampVec = linspace(currVec(1),currVec(3),floor((currVec(3) - currVec(1))/currVec(2))+1); 
	%Normalize
	nratios = ratios/max(abs(ratios)); 		
	ntime = time/timeRes; 
	pulseTime = sum(ntime); %number of 50 us segments in single waveform for use in PLI
	%Make whole number version of PL
	plVecInt = [];
	for i=ampVec
		for j=1:numel(nratios)
			for k=1:ntime(j)
				plVecInt = [plVecInt nratios(j)*i];
			end
		end
	end
	%Make PL
	plVec = daqify(plVecInt, range,numBits);
	plVec = [repmat(preDaqVec, size(plVec,1), 1) plVec];
	plVec = [setRecordVec; plVec];
	pl = int16(bi2de(plVec));

	%Make PLI
	pli = [1:extraCommandNum linspace(0,numel(plVecInt),numel(ampVec)+1)+1+extraCommandNum]; %assume each command is only 50 us
	pli = int32(pli);

	%Make ESF
	chanVec = 1:numChannel;
	esfmat = zeros(iter*numel(ampVec)*numChannel + extraCommandNum*numChannel, 3);
	%Set all channels to record
	for i=chanVec
		esfmat(i,:) = [0, i, 1];
	end
	%random scan
	if randflg == 1
		esfmatInd = numChannel + 1;
		timInc = timeRes;
		for i=1:length(pli) %same as numel(ampVec)
			for j=1:iter
				nchanVec = chanVec(randperm(length(chanVec)));
				for k=1:numChannel
					esfmat(esfmatInd,:) = [timInc, nchanVec(k), i];
					esfmatInd = esfmatInd + 1;
					timInc = timInc + timeRes + gapParam;
				end
			end
		end
	end
	esfmat = int32(esfmat);
	
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
