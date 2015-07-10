function fitGMLM = scramble_fitGMLM(fitGMLM)
nFrontEnds = length(fitGMLM.Linear.filter);

if(nFrontEnds ==1)
return 
end

% Random orthogonal transform..
% Find filter matrix
sdim = length(fitGMLM.Linear.filter{1});
Filters=zeros(nFrontEnds,sdim);

for ifilt=1:nFrontEnds
Filters(ifilt,:) = fitGMLM.Linear.filter{ifilt};    
end

[q,r] = qr(randn(nFrontEnds,nFrontEnds));
Filtersnew = q*Filters; 

for ifilt=1:nFrontEnds
fitGMLM.Linear.filter{ifilt} = Filtersnew(ifilt,:)';     
end
end