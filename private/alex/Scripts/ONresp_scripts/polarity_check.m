clear
load('/Users/alexth/Desktop/old_stuff/ONresp/datarun.mat')

tmp=zeros(500,4,size(datarun,2));


for i=1:size(datarun,2)
    k=[datarun(i).lf.nd7 datarun(i).lf.control];
    k=k-repmat(mean(k),500,1);
    k=k./abs(repmat(max(k)-min(k),500,1));
    tmp(:,:,i)=k;
end

tmp(isnan(tmp))=0;
plot(k)

for i=1:4
    t=squeeze(tmp(50:250,i,:));
    x=t';
    x=cov(x);
    [V,~]=eig(x);
    pc_vectors=V(:,end);
    pc1(:,i)=t'*pc_vectors;
end


onCells=find(sum(pc1<0,2)>2); % ON cells
offCells=find(sum(pc1>0,2)>2); % OFF cells

strangeCells=find(sum(pc1>0,2)==2); % changing cells? Check:
cnt=1;
for i=strangeCells'
    subplot(3,4,cnt)
    cnt=cnt+1;
    plot(tmp(:,:,i))
end
clear tmp
tmp(onCells)=1;
tmp(offCells)=-1;
tmp(strangeCells)=0;

for i=1:size(datarun,2)
    datarun(i).polarity=tmp(i);
end

save('/Users/alexth/Desktop/old_stuff/ONresp/datarun.mat','datarun')


load('/Users/alexth/Desktop/old_stuff/ONresp/datarun_nd.mat')
for i=1:size(datarun_nd,2)
    datarun_nd(i).polarity=tmp(i);
end
save('/Users/alexth/Desktop/old_stuff/ONresp/datarun_nd.mat','datarun_nd')
