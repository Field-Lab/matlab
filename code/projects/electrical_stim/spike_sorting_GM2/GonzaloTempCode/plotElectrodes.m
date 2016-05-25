function [aa elcenter dist distStim]=plotElectrodes(templates,thresEI,varargin)
load arrayPositions512

%positions=(positions-repmat(mean(positions),size(positions,1),1))./repmat(std(positions),size(positions,1),1)

els=[];
for n=1:length(templates)
  
    [a b]=sort(max(abs(templates{n}')),'descend');
  elcenter(n,:)=nanmean(positions(b(1:3),:));
end

for n1=1:length(templates)
    for n2=1:length(templates)
        dist(n1,n2)=norm(elcenter(n1,:)-elcenter(n2,:));
    end
end

if(nargin==3)
    stimElec=varargin{1};
    ind=setdiff([1:512],stimElec);
    for n=1:length(templates)
        distStim(n)=norm(positions(stimElec,:)-elcenter(n,:));
    end
else
    ind=[1:512];
end


els=[];
for n=1:length(templates)
    templates{n}=templates{n}(ind,:);
    [a b]=sort(max(abs(templates{n}')),'descend');
    ind2=find(a>thresEI);
    if(isempty(ind2))
       ind2=1;
    end
   elsLocal{n}=b(ind2);
   sizesLocal{n}=a(ind2);
   vec(n,:)=zeros(1,length(ind));
   vec(n,elsLocal{n})=sizesLocal{n};

     els=union(b(ind2),els);
end

figure(1000)
subplot(2,2,1)
scatter(elcenter(:,2),elcenter(:,1));
hold on
for n=1:length(templates)
    text(elcenter(n,2)+0.025,elcenter(n,1)+0.05,num2str(n))
end
aa=vec*vec';
subplot(2,2,4)
imagesc(aa)
subplot(2,2,3)
imagesc(dist)

