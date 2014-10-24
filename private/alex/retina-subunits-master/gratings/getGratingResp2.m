function g = getGratingResp(grat,celldat)

%-------------------------------------------
%
% g = getGratingResp(grat,dat,fit)
%
% extract some useful response parameters from
% a grating run, for a single cell. alternatively,
% if 'fit' is provided, extract the same parameters
% for a model cell
%
% grat should contain the fields:
% grat.spikes -- cell array of spike times
% grat.stimulus -- structure with stimulus info
% grat.cell_ids -- vector of cell ids
% grat.cone_inputs -- matrix of cone generator signals
%                     (cones x time points x conditions)
%
% dat should contain the fields:
% dat.coneIds -- vector of cone numbers for the cell, cones x 1
% dat.rgcId -- what cell are we doing the analysis on?
%
% fit is an optional structure, 
% see 'fitA.m' and 'fitB.m' for more info
%
% freeman, 5-5-2012
%-------------------------------------------

% set some constants
REFRESH = 60; % was 120 for 2010-03-05-2/data016
STIMDUR = 8;

% get spike times for this cell
if sum(grat.cell_ids==celldat.rgcId)
    spikeTimes = grat.spikes{grat.cell_ids==celldat.rgcId};
else
    fprintf('(getGratingResp) no matching cell for %g\n',celldat.rgcId);
    g = [];
    return
end

% set some timing parameters
bin = 1/REFRESH;
nRep = grat.stimulus.repetitions * 2;
cycleDur = grat.stimulus.params.TEMPORAL_PERIOD;
nCycle = (STIMDUR*REFRESH)/cycleDur;

% get number of periods
periods = sort(unique([grat.stimulus.combinations.SPATIAL_PERIOD]),'descend');
perInds = find([grat.stimulus.combinations.SPATIAL_PERIOD]==periods(1));
% get spatial frequencies
perVals = (1000/5.8)./fliplr(grat.stimulus.params.SPATIAL_PERIOD);
% get number of unique phases for one period (ASSUMING SAME FOR ALL)
phases = unique([grat.stimulus.combinations(perInds).SPATIAL_PHASE]);
nPer = length(periods);
nPh = length(phases);
% preallocate matrices
coneMaxFrame = zeros(size(celldat.coneIds,1),nPer,nPh);
timeCourses = zeros(cycleDur*nCycle,nRep,nPer,nPh);
stim = zeros(cycleDur,nPer,nPh,length(celldat.coneIds));
cycleAve = zeros(cycleDur,nPer,nPh);
f1 = zeros(nRep,nPer,nPh);
f2 = zeros(nRep,nPer,nPh);

% get the time course for each reptition, spatial period, and phase
for iper=1:nPer
    perInds = find([grat.stimulus.combinations.SPATIAL_PERIOD]==periods(iper));
    phases = unique([grat.stimulus.combinations(perInds).SPATIAL_PHASE]);
    for iph=1:length(phases)
        trialInd = find([grat.stimulus.combinations.SPATIAL_PERIOD]==periods(iper) & ...
            [grat.stimulus.combinations.SPATIAL_PHASE]==phases(iph));
        coneMaxFrame(:,iper,iph) = grat.cone_inputs(celldat.coneIds,8,trialInd(1));
        repInds = [find(grat.stimulus.trial_list == trialInd(1)), find(grat.stimulus.trial_list == trialInd(2))];
        for irep=1:length(repInds)
            inds = linspace(0,STIMDUR,(STIMDUR/bin)+2);
            t = hist(spikeTimes,inds+grat.stimulus.triggers(repInds(irep)));
            t = t(2:end-1);   
            tmpStim = grat.cone_inputs(celldat.coneIds,:,trialInd);
            if irep == 1
                stim(:,iper,iph,:) = tmpStim(:,1:cycleDur)'; 
            end
            timeCourses(:,irep,iper,iph) = t;            
        end
    end
end

% get the cycle average response
for iper=1:nPer
    for iph=1:nPh
        tmp = reshape(timeCourses(:,:,iper,iph),[cycleDur nCycle*nRep]);
        cycleAve(:,iper,iph) = mean(tmp,2);
    end
end

% get the f1 and f2 for each repitition, period, and phase
for iper=1:nPer
    for iph=1:nPh
        for irep=1:3
            rr = timeCourses(:,irep,iper,iph);
            nfft=2^nextpow2(length(rr));
            y=fft(rr,nfft)/length(rr);
            f=1/bin/2*linspace(0,1,nfft/2);
            [~, t1]=min(abs(f-REFRESH/cycleDur));
            [~, t2]=min(abs(f-REFRESH/cycleDur*2));
            f1(irep,iper,iph)=2*abs(y(t1));
            f2(irep,iper,iph)=2*abs(y(t2));
        end
    end
end

timeCourseNorm = celldat.timeCourse / norm(celldat.timeCourse);
timeCourseNorm = flipud(resample(flipud(timeCourseNorm),8,1));
timeCourseNorm = timeCourseNorm / norm(timeCourseNorm);
for iPer=1:nPer
    for iPh=1:nPh
        stim(:,iPer,iPh,:) = tempFilterCones(squeeze(stim(:,iPer,iPh,:)),timeCourseNorm,'circular');
    end
end

% collect results
g.timeCourses = timeCourses;
g.cycleAve = cycleAve;
g.stim = stim;
g.f1 = f1;
g.f2 = f2;
g.perVals = perVals;
g.coneMaxFrame = coneMaxFrame;
g.nPer = nPer;
g.nPh = nPh;
g.nRep = nRep;
g.cycleDur = cycleDur;
g.nCycle = nCycle;