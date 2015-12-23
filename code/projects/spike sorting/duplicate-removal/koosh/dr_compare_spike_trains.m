function [merged] = dr_compare_spike_trains(s, e, duration, threshold, samplingRate)

% cross correlate the two spike trains
binned1=histc(s.a/samplingRate,0:10/1000:duration-10/1000);
binned2=histc(s.b/samplingRate,0:10/1000:duration-10/1000);

r=corrcoef(binned1,binned2);
r=r(1,2);

% if the clusters are on the same electrode we use a smaller search
% window to see if spikes are "identical" than if they are not
if e.a == e.b
    multiple = 1;
else
    multiple = 4;
end

% 5 is the minimum number of samples between spikes allowed by vision
% our threshold should be some multiple of this number (.25 ms/5 samples)
theta = 5*multiple;

if max(s.a) > max(s.b)
    temp = s.a; eTemp = e.a;
    s.a = s.b; e.a = e.b;
    s.b = temp; e.b = eTemp;
end
    
% if they are "identical" then combine the spikes
if r >=threshold

    % kill duplicate spikes
    bPos = 1;
    for aPos=1:length(s.a)
        % check if we are beyond the beginning of the search window
        while s.b(bPos) > s.a(aPos)-theta
            if bPos == 1
                break;
            end
            bPos = bPos - 1; 
        end
            
        % move iterator to beginning of search window if necessary
        % since max(s.b) must be greater than max(s.a) this location must
        % be before the end of spike train b
        while s.b(bPos) < s.a(aPos)-theta
            bPos = bPos + 1;
        end
        
        while true
            if bPos > length(s.b) || s.b(bPos) > s.a(aPos)+theta
                break;
            end
            
            % kill all bad samples in the search window
            s.b(bPos) = [];
            % don't increment the pointer because the vector 
            % is now shorter by one sample
        end
        
        if bPos > length(s.b)
            break;
        end
    end
    
    merged.s = unique(union(s.a,s.b));
    merged.e = e.a;
else
    merged = [];
end