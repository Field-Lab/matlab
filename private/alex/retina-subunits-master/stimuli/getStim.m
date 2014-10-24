function [train test] = getStim(simParams,dat,train,test)

if isempty(simParams.len) && exist('train','var')
    simParams.len(1) = size(train.X_ct,2);
	simParams.len(2) = size(test.X_ct,2);
end

if ~isfield(simParams,'C')
    simParams.C = size(train.X_ct,1);
end

dimsTrain = [simParams.C simParams.len(1)];
dimsTest = [simParams.C simParams.len(2)];

if strcmp(simParams.type,'sinusoidal') || strcmp(simParams.type,'jitter') || ...
    strcmp(simParams.type,'smooth') || strcmp(simParams.type,'smoothThresh') || ...
    strcmp(simParams.type,'gaussianCones')
        
    mnValX = min(dat.locationsFit(:,1));
    mxValX = max(dat.locationsFit(:,1));
    mnValY = min(dat.locationsFit(:,2));
    mxValY = max(dat.locationsFit(:,2));
    res = 100;
    xRange = linspace(mnValX-1,mxValX+1,res);
    yRange = linspace(mnValY-1,mxValY+1,res);
    coneVar = 0.05;
    [X Y] = meshgrid(xRange,yRange);
    [distance diff] = getDistances(dat.locationsFit,[X(:) Y(:)]);
    diff = reshape(permute(diff,[1 3 2]),[dat.nConesFit*res^2,2]);
    cones = evalGaussian(diff,[1 0; 0 1]*coneVar);
    cones = reshape(cones,[dat.nConesFit res^2]);
    cones = normMat(cones,'L2',2);
    M = cones;
    
    train.X_ct = zeros(dimsTrain);
    test.X_ct = zeros(dimsTest);
end

switch simParams.type
    case 'uniform'
        train.X_ct = rand(dimsTrain)*2-1;
        test.X_ct = rand(dimsTest)*2-1;
    case 'gaussian'
        train.X_ct = randn(dimsTrain);
        test.X_ct = randn(dimsTest);
    case 'binary'
        train.X_ct = real(rand(dimsTrain)>0.5)*2-1;
        test.X_ct = real(rand(dimsTest)>0.5)*2-1;
    case 'original'
        train.X_ct = train.X_ct(1:dimsTrain(1),:);
        test.X_ct = test.X_ct(1:dimsTest(1),:);
    case 'correlated'
        [dist diff] = getDistances(dat.locationsFit,dat.locationsFit);
        diff = reshape(permute(diff,[1 3 2]),[dat.nConesFit*dat.nConesFit,2]);
        gaussVals = evalGaussian(diff,[1 0; 0 1]*simParams.sd,[],1);
        gaussVals = reshape(gaussVals,[dat.nConesFit dat.nConesFit]);
        train.X_ct = randn(dimsTrain)*sqrtm(gaussVals);
        test.X_ct = randn(dimsTest)*sqrtm(gaussVals);
    case 'sinusoidal'
        
        % cycle over different frequencies

        switch simParams.grating
            case 'drifting'
                blockLen = 16;
            case 'modulated'
                blockLen = 32;
        end
        
        nFreqs = length(simParams.freq);
        
        nBlocksTrain = floor(dimsTrain(2)/blockLen);
        blocksPerFreqTrain = floor(nBlocksTrain/nFreqs);
        tmp = repmat(simParams.freq,[blockLen*blocksPerFreqTrain 1]);
        freqsTrain = tmp(:);
        freqsTrain(end+1:dimsTrain(2)) = freqsTrain(end);
        
        nBlocksTest = floor(dimsTest(2)/blockLen);
        blocksPerFreqTest = floor(nBlocksTest/nFreqs);
        tmp = repmat(simParams.freq,[blockLen*blocksPerFreqTest 1]);
        freqsTest = tmp(:);
        freqsTest(end+1:dimsTest(2)) = freqsTest(end);
        
        switch simParams.grating
            case 'drifting'
                thisAmp = 1;
                phases = linspace(0,2*pi,16);
                phases = phases(1:end-1);
            case 'modulated'
                thisPhase = simParams.phase;
                amps = linspace(-1,1,16);
                amps = [amps fliplr(amps(2:end-1))];
        end
        for i=1:dimsTrain(2)
            switch simParams.grating
                case 'drifting'
                    thisPhase = phases(mod(i,length(phases))+1);
                case 'modulated'
                    thisAmp = amps(mod(i,length(amps))+1);
            end
            s = thisAmp*cos((2*pi)*freqsTrain(i)*X+(2*pi)*0*Y+thisPhase);
            c = M*s(:);
            train.X_ct(:,i) = c;
        end
        for i=1:dimsTest(2)
            switch simParams.grating
                case 'drifting'
                    thisPhase = phases(mod(i,length(phases))+1);
                case 'modulated'
                    thisAmp = amps(mod(i,length(amps))+1);
            end
            s = thisAmp*cos((2*pi)*freqsTest(i)*X+(2*pi)*0*Y+thisPhase);
            c = M*s(:);
            test.X_ct(:,i) = c;
        end
    case {'smooth','smoothThresh'}
        
        for i=1:dimsTrain(1)
            if length(simParams.f) == 1
                s = mkFract(res,simParams.f);
            else
                tmpF = shuffle(simParams.f);
                tmpF = tmpF(1);
                s = mkFract(res,tmpF);
            end
            c = s(:)'*M;
            train.X_ct(i,:) = c;
        end
        for i=1:dimsTest(1)
            if length(simParams.f) == 1
                s = mkFract(res,simParams.f);
            else
                tmpF = shuffle(simParams.f);
                tmpF = tmpF(1);
                s = mkFract(res,tmpF);
            end
            c = s(:)'*M;
            test.X_ct(i,:) = c;
        end
        if strcmp(simParams.type,'smoothThresh')
            tmp = zeros(size(train.X_ct));
            tmp(train.X_ct<-2) = -1;
            tmp(train.X_ct>2) = 1;
            train.X_ct = tmp;
            tmp = zeros(size(test.X_ct));
            tmp(test.X_ct<-2) = -1;
            tmp(test.X_ct>2) = 1;
            test.X_ct = tmp;
        end 
        %train.X_ct = train.X_ct/std(train.X_ct(:));
        %test.X_ct = test.X_ct/std(test.X_ct(:));
    case 'jitter'
        stixSz = simParams.stixSz; % stixel size
        for i=1:dimsTrain(1)
            s = randn(round(res/stixSz));
            s = imresize(s,[res res],'nearest');
            shiftX = ceil(rand*stixSz-1);
            shiftY = ceil(rand*stixSz-1);
            s = circshift(s,[shiftX shiftY]);
            c = s(:)'*M;
            train.X_ct(i,:) = c;
        end
         for i=1:dimsTest(1)
            s = randn(round(res/stixSz));
            s = imresize(s,[res res],'nearest');
            shiftX = ceil(rand*stixSz-1);
            shiftY = ceil(rand*stixSz-1);
            s = circshift(s,[shiftX shiftY]);
            c = s(:)'*M;
            test.X_ct(i,:) = c;
         end
         
    case 'gaussianCones'
        train.X_ct = randn(res^2,dimsTrain(1))'*M;
        test.X_ct = randn(res^2,dimsTest(2))'*M;
end

if simParams.p < 1
train.X_ct(rand(size(train.X_ct))>simParams.p) = 0;
test.X_ct(rand(size(test.X_ct))>simParams.p) = 0;
end