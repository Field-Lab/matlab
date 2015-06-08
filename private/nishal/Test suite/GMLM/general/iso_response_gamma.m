function iso_response_gamma(binnedResponsesbigd,maskedMov,fitGMLM,gamma,su)
%% Extract filters
filters = fitGMLM.Linear.filter;
nFrontEnds = length(su);
mu=fitGMLM.mu;


mov_filtered=maskedMov;
%% Calculate sub-unit response
kx=cell(nFrontEnds,1);
kkx = cell(nFrontEnds,1);

for ifilter=1:2%nFrontEnds
    kkx{ifilter}= filters{su(ifilter)}'*mov_filtered;
    kx{ifilter} = (kkx{ifilter}.*(kkx{ifilter}>0)).^gamma;
end

%% Empirical sub-unit activation distribution
dimensions=[]; binEdge_log=cell(2,1);
[f,xi]=ecdf(kx{1}); eps=0.05
binEdge=[];
for it=0:eps:1
binEdge = [binEdge;xi(find(f<=it,1,'last'))];
end
binEdge=unique(binEdge);
binEdge_log{1}=binEdge;

[f,xi]=ecdf(kx{2}); eps=0.05
binEdge=[];
for it=0:eps:1
binEdge = [binEdge;xi(find(f<=it,1,'last'))];
end
binEdge=unique(binEdge);
binEdge_log{2}=binEdge; 
Ymesh = repmat(binEdge_log{1},[1,length(binEdge_log{2})]);

Xmesh = repmat(binEdge_log{2}',[length(binEdge_log{1}),1]);
for isu=1:nFrontEnds
    
su_range{isu} =binEdge_log{isu};%([(min(kx{isu})):3:(max(kx{isu}))]); %Use freedman - diaconis binning .. new matlab has it .. 
dimensions = [dimensions,length(su_range{isu})];
end

grid=repmat(0,dimensions);
grid_total=repmat(0,dimensions);

for itime=1:length(binnedResponsesbigd)
   
    su_id = cell(nFrontEnds,1);
    for isu=1:nFrontEnds
        idx=1:length(su_range{isu});
        su_id{isu} = idx(find(kx{isu}(itime)>=su_range{isu},1,'last'));
    end
    grid(su_id{1},su_id{2})= grid(su_id{1},su_id{2}) + binnedResponsesbigd(itime);
   
    grid_total(su_id{1},su_id{2}) =  grid_total(su_id{1},su_id{2}) +1;
end
    
grid=grid./grid_total;

%% Plot
  figure;
  contourf(Xmesh,Ymesh,grid_total);
  figure; 
  contourf(Xmesh,Ymesh,grid,50);
end