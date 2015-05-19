function iso_response_bias_gmlm(binnedResponsesbigd,maskedMov,fitGMLM)
%% Extract filters
filters = fitGMLM.Linear.filter;
bias = fitGMLM.Linear.bias;
nFrontEnds = length(filters);
mu=fitGMLM.mu;


mov_filtered=maskedMov;
%% Calculate sub-unit response
kx=cell(nFrontEnds,1);

for ifilter=1:nFrontEnds
    kx{ifilter} = exp(filters{ifilter}'*mov_filtered + bias{ifilter});
end


%% Empirical sub-unit activation distribution
dimensions=[];
[f,xi]=ecdf(kx{1}); eps=0.05;
binEdge=[];
for it=0:eps:1
binEdge = [binEdge;xi(find(f<=it,1,'last'))];
end
Xmesh = repmat(binEdge',[length(binEdge),1]);
Ymesh = repmat(binEdge,[1,length(binEdge)]);

for isu=1:nFrontEnds
    
su_range{isu} =binEdge;%([(min(kx{isu})):3:(max(kx{isu}))]); %Use freedman - diaconis binning .. new matlab has it .. 
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
  contourf(Xmesh,Ymesh,grid_total,50);
  figure; 
   contourf(Xmesh,Ymesh,grid,50);
end