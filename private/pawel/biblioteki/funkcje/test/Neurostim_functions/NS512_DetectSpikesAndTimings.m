function Events=NS512_DetectSpikes(DataTraces,ThresholdNegative,ThresholdPositive);

mins=min(DataTraces,[],3);
EventsMin=floor((-sign(mins+ThresholdNegative)+1)/2); %2-D array - which traces for each electrode go below -ThresholdNegative
maxs=max(DataTraces,[],3);
EventsMax=floor((sign(maxs-ThresholdPositive)+1)/2); %2-D array - which traces for each electrode go above ThresholdPositive
find(EventsMax==1);

Events=EventsMin.*EventsMax; %2-D array - combination of the two arrays above