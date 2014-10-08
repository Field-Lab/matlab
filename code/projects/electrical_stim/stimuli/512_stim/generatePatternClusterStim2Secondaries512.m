function [arrayForCluster smallAmpPatterns] = generatePatternClusterStim2Secondaries512(clusterElectrodes, centerElec, relAmps)

% generates array of patterns using pairs (one primary, one secondary electrode) and triplets (one
% primary, 2 secondary electrodes) of electrodes with all secondary electrode polarity combinations and relative
% amplitude combinations.  An additional set of relative amplitudes can be chosen to apply only in
% electrode pairs (1 primary, 1 secondary)
%
% the current output of the primary electrode is always set to 1, and all secondary amplitudes are
% specified relative to the primary amplitude.  Actual current amplitudes are specified in the
% labview software.
%
% arguments:
%             clusterElectrodes - vector of electrode numbers (all electrodes in a single cluster)
%                    centerElec - electrode number corresponding to primary electrode
%                relAmps.normal - secondary electrode amplitudes, relative to the primary electrode
%                                 amplitude, to be used in all combinations of 2 and 3 electrodes
%             relAmps.pairsOnly - secondary electrode amplitudes, relative to the primary
%                                 electrode amplitude, to be used in combinations of 2 electrodes only (primary + 1 secondary)
%        relAmps.pairsOnlySmall - same as .pairsOnly, except that these patterns are played in a separate movie chunk from the rest to avoid
%                                 rounding issues that occur when very low and very high currents are used on the same electrode in the same movie chunk%                                 
%
% outputs:
%               arrayForCluster - array of values representing relative amplitudes of pulses on each
%                                 electrode, with first dimension corresponding to electrodes and 
%                                 second dimension corresponding to different patterns
%              smallAmpPatterns - binary vector flagging patterns that only
%                                 contain relative secondary amplitudes smaller than the value
%                                 'smallAmpThresh'

% author: Lauren Hruby Jepson
% last checked 2010-10-13
% added pseudo-local return stim patterns and anodal primary-alone stimulus 2010-11-10
% added smallAmpThresh/smallAmpPatterns 2011-02-28
%
% 2011-04-01 removied smallAmpThresh, replaced with relAmps.pairsOnlySmall,
% which specifies directly which amplitudes to separate out
%
% 2011-06-21 changed amplitude for secondary-alone stimuli to be equal to
% the maximum relative amplitude used (typically, 1.5 for 30 micron and 2
% for 60 micron)


nElec = length(clusterElectrodes);
nRelAmps = length(relAmps.normal);
surrElecInd = find(clusterElectrodes ~= centerElec);
centerElecInd = find(clusterElectrodes == centerElec);

% determines which of the secondary electrodes are next to eachother and returns a list of such
% pairs (should be <= pairs)
combinations = nearestNeighborsFinder512(clusterElectrodes(surrElecInd));

arrayForCluster = [];

relAmpsAll = [relAmps.pairsOnly relAmps.pairsOnlySmall relAmps.normal];
maxRelAmp = max(relAmpsAll);

%% pairs patterns
for j = 1:length(relAmpsAll)
    for k = 1:nElec-1
        %both +
        arrayForCluster = [arrayForCluster zeros(nElec,1)];%#ok<AGROW> %add column of zeros to array
        arrayForCluster(surrElecInd(k), end) = relAmpsAll(j);
        arrayForCluster(centerElecInd, end) = 1;

        %primary +, secondary -
        arrayForCluster = [arrayForCluster zeros(nElec,1)];%#ok<AGROW> %add column of zeros to array
        arrayForCluster(surrElecInd(k), end) = -1*relAmpsAll(j);
        arrayForCluster(centerElecInd, end) = 1;
    end
end

%% triplets patterns

%generates list of all amplitude combinations on pairs of secondary electrodes
ampCombs = zeros(nRelAmps^2, 2);
for i = 1:nRelAmps
    for j = 1:nRelAmps
        ampCombs((i-1)*nRelAmps+j,1) = relAmps.normal(i);
        ampCombs((i-1)*nRelAmps+j,2) = relAmps.normal(j);
    end
end

for j = 1:size(ampCombs,1)
    for k = 1:size(combinations,1)

        %all +
        arrayForCluster = [arrayForCluster zeros(nElec,1)];%#ok<AGROW> %add column of zeros to array
        arrayForCluster(surrElecInd(combinations(k,1)), end) = ampCombs(j,1);
        arrayForCluster(surrElecInd(combinations(k,2)), end) = ampCombs(j,2);
        arrayForCluster(centerElecInd, end)                  = 1;

        %p+, s1+, s2-
        arrayForCluster = [arrayForCluster zeros(nElec,1)];%#ok<AGROW> %add column of zeros to array
        arrayForCluster(surrElecInd(combinations(k,1)), end) = ampCombs(j,1);
        arrayForCluster(surrElecInd(combinations(k,2)), end) = -1*ampCombs(j,2);
        arrayForCluster(centerElecInd, end)                  = 1;

        %p+, s1-, s2+
        arrayForCluster = [arrayForCluster zeros(nElec,1)];%#ok<AGROW> %add column of zeros to array
        arrayForCluster(surrElecInd(combinations(k,1)), end) = -1*ampCombs(j,1);
        arrayForCluster(surrElecInd(combinations(k,2)), end) = ampCombs(j,2);
        arrayForCluster(centerElecInd, end)                  = 1;

        %p+, s1-, s2-
        arrayForCluster = [arrayForCluster zeros(nElec,1)];%#ok<AGROW> %add column of zeros to array
        arrayForCluster(surrElecInd(combinations(k,1)), end) = -1*ampCombs(j,1);
        arrayForCluster(surrElecInd(combinations(k,2)), end) = -1*ampCombs(j,2);
        arrayForCluster(centerElecInd, end)                  = 1;
    end 
end


%% imitation local return stimulus
arrayForCluster = [arrayForCluster zeros(nElec,1)];%#ok<AGROW> %add column of zeros to array
arrayForCluster(surrElecInd, end) = -1/6;
arrayForCluster(centerElecInd, end) = 1;

%% single electrode stimuli

% primary electrode alone (both polarities)
arrayForCluster = [arrayForCluster zeros(nElec,1)];%add column of zeros to array
arrayForCluster(centerElecInd, end) = 1;

arrayForCluster = [arrayForCluster zeros(nElec,1)];%add column of zeros to array
arrayForCluster(centerElecInd, end) = -1;

% secondary electrodes alone (both polarities)
for i = 1:nElec-1
    arrayForCluster = [arrayForCluster zeros(nElec,1)];%#ok<AGROW> %add column of zeros to array
    arrayForCluster(surrElecInd(i), end) = maxRelAmp;
    
    arrayForCluster = [arrayForCluster zeros(nElec,1)];%#ok<AGROW> %add column of zeros to array
    arrayForCluster(surrElecInd(i), end) = -1*maxRelAmp;
end

%% flag all patterns that contain relative amplitudes equal to what's
% listed in relAmps.pairsOnlySmall

smallAmpPatterns = false(1,size(arrayForCluster,2));

for ii = 1:size(arrayForCluster, 2)
    sElecAmps = arrayForCluster(surrElecInd,ii);
    for jj = 1:length(relAmps.pairsOnlySmall)
        if any(abs(sElecAmps) == relAmps.pairsOnlySmall(jj))
            smallAmpPatterns(ii) = true;
            break
        end
    end
end



% for ii = 1:size(arrayForCluster,2)
%     if any((arrayForCluster(:,ii)~=0 & abs(arrayForCluster(:,ii))<smallAmpThresh))
%         sElecAmps = arrayForCluster(surrElecInd,ii);
%         %remove zero amplitudes
%         sElecAmps = sElecAmps(sElecAmps ~= 0);
%         if all(abs(sElecAmps)<smallAmpThresh) %only flag if all nonzero secondary amplitudes are below threshold
%             smallAmpPatterns(ii) = true;
%         end
%     end
% end




