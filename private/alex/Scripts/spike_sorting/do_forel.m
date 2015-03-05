function do_forel

global thr rgc uicontr_units spikes getXdata getYdata visibleUnits signUnit fatPoints ifUpdateISI

% select subset of spikes and align by minima
for j=1:size(rgc,2)
    if get(uicontr_units(j),'Value')==1
        subset=j;
        break
    end
end

spikesSelected=spikes(:,rgc{subset});
spikesSubset=zeros(1,1:size(spikes,2));
spikesSubset(rgc{subset})=1;
% spikesSubset=min(spikes)<thr;
% spikesSelected=spikes(:,spikesSubset);


[~,pos]=min(spikesSelected(10:14,:));
spikes_tmp=zeros(20,size(spikesSelected,2));
pos=pos+9;
for i=10:14
    first_point=i-9;
    last_point=i+10;
    spikes_tmp(:,pos==i)=spikesSelected(first_point:last_point,pos==i);
end

clear pos first_point last_point spikesSelected

% do forel
sparseness=ceil(size(spikes_tmp,2)/2000);
NofClusters=10;
cluster=cell(sparseness,NofClusters);
realTime=tic;
for share=1:sparseness
    if share==1
        estimTime=tic;
    elseif share==2
        estimTime=toc(estimTime);
        fprintf('\nEstimated waiting time: %d s', ceil(estimTime*sparseness))
    end

    m=spikes_tmp(:,share:sparseness:end)';
    mm=pdist(m);
    
    common=zeros(size(m,1));
    common = squareform(mm);
    %     clear mm m
    common=common.*common;
    spread=zeros(size(common,1));
    tmp=std(common);
    spread=repmat(tmp,size(tmp,2),1);
    common=common-spread;
    common(common>0)=0;
    common(common<0)=1;
    
    k=1;
    while k<NofClusters+1
        
        if ~isempty(cluster{share,1})
            common(central_element(m),:)=0;
            common(:,central_element(m))=0;
        end
        if nnz(common)==0
            break
        end
        posElem=sum(common>0);
        [~,b]=max(posElem);
        central_element=find(common(:,b));
        commonTMP=common(central_element,central_element);
        f=sum(commonTMP);
        m=find(f>mean(f));
        
        cluster(share,k)={central_element(m)};
        k=k+1;
    end    
end
realTime=toc(realTime);
fprintf('\nreal time: %d s\n',ceil(realTime))

% make timing independent of sparseness
all_spike_times=1:size(spikes_tmp,2);
for share=1:sparseness
    k=all_spike_times(share:sparseness:end);
    for i=1:NofClusters
        cl_tmp{i,share}=k(cluster{share,i});
    end
end


%%%%%%%%%%%%%%%%%% GETTING CORRELATED CLUSTERS %%%%%%%%%%%%%%%%%%%

templates=zeros(NofClusters*sparseness,20); % waveform lenth by amount of clusters per sorting by number of sortings
cnt=1;
for share=1:sparseness
    for i=1:NofClusters
        templates(cnt,:)=mean(spikes_tmp(:,cl_tmp{i,share}),2);
        cnt=cnt+1;
    end
end

clear ordering
common=pdist(templates);
common = squareform(common);
m=max(common(:));
i=1;cnt=1;kk=[];
lisst=1:NofClusters*sparseness;
while i<NofClusters*sparseness
    if sum(diff(common(i,:)))~=0
        a=find(common(i,:)<200);
        if length(a)>=sparseness/2
            kk=sort([kk a]);
            lisst(2,a)=cnt;
            ordering{cnt}=a;
            common(a,:)=m;
            common(:,a)=m;
            cnt=cnt+1;
        end
    end
    i=i+1;
end


for i=1:length(ordering)
    clusters{i}=sort([cl_tmp{ordering{i}}]);
end

% assigning to rgc
tmp=find(spikesSubset);
for i=1:min(length(clusters),9)
    clusters{i}=tmp(clusters{i});
end


rgc=clusters;

visibleUnits=ones(1,size(rgc,2));
signUnit=zeros(1,size(rgc,2));
fatPoints=zeros(1,size(rgc,2));
ifUpdateISI=1;
redraw(get(getXdata,'value'),get(getYdata,'value'))