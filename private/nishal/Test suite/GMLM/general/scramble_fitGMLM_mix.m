function fitGMLM = scramble_fitGMLM_mix(fitGMLM)
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

frac = rand(1);
[q,r] = qr(randn(nFrontEnds,nFrontEnds));
e=ones(nFrontEnds,1); e=e/norm(e);
q2 = (eye(nFrontEnds) - e*e')*q+ e*e';
Filtersnew = q2*Filters; 

for ifilt=1:nFrontEnds
fitGMLM.Linear.filter{ifilt} = Filtersnew(ifilt,:)';     
end
end