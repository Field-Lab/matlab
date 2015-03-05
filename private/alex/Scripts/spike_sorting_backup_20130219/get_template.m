function get_template(which_units)

global spikes rgc distance2template

% align by min
[~,pos]=min(spikes(10:14,:));
spikes_tmp=zeros(20,size(spikes,2));
pos=pos+9;
for i=10:14
    first_point=i-9;
    last_point=i+10;
    spikes_tmp(:,pos==i)=spikes(first_point:last_point,pos==i);
end

% normalize amplitude (from 0 to -1)
spikes_tmp=spikes_tmp-repmat(max(spikes_tmp),20,1);
spikes_tmp=-spikes_tmp./repmat(min(spikes_tmp),20,1);

clear pos

template=zeros(20,length(which_units));

distance2template=zeros(size(spikes,2),size(template,2));

%get templates
cnt=1;
for k=which_units    
    tmp=spikes_tmp(:,rgc{k});
    
%     % normalize amplitude (from 0 to -1)
%     tmp=tmp-repmat(max(tmp),20,1);
%     tmp=-tmp./repmat(min(tmp),20,1);    
    
    sparseness=ceil(size(tmp,2)/5000); % take 5000 equally spaced spikes from unit
    tmp=tmp(:,1:sparseness:end);
    
    common=pdist(tmp');
    common = squareform(common);
    [common,IX]=sort(common);
    [~,a]=min(mean(common(1:min(size(tmp,2)/2,2500),:))); % take 2500 closest spikes from the most central spike
    a=IX(:,a);
    template(:,cnt)=mean(tmp(:,a(1:min(size(tmp,2)/2,2500))),2);
    
    tt=repmat(template(:,cnt),1,size(spikes_tmp,2));
    distance2template(:,cnt)=sqrt(sum((tt-spikes_tmp).*(tt-spikes_tmp)));
    
    cnt=cnt+1;
end

col=[1,0,0;0,1,0;0.5,0.5,0;0.25,0.75,0.75;0.75,0.25,0.75;0.05,0.58,0.05;0.9,0.45,0.1;0.4,0.25,0.75;0.86,0.86,0];
figure
for i=1:size(template,2)
    plot(template(:,i),'color',col(which_units(i),:),'linewidth',3);
    hold on
end