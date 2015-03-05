function subset_pca(checked)

global spikes rgc isi_all timeUnit pc1 tdata

spikes_tmp=spikes(:,rgc{checked});
tvector=tdata(rgc{checked});
% align by min
[~,pos]=min(spikes_tmp(10:14,:));
norm_spikes=zeros(20,size(spikes_tmp,2));
pos=pos+9;
for i=10:14
    first_point=i-9;
    last_point=i+10;
    norm_spikes(:,pos==i)=spikes_tmp(first_point:last_point,pos==i);
end

% normalize amplitude (from 0 to -1)
norm_spikes=norm_spikes-repmat(max(norm_spikes),20,1);
norm_spikes=-norm_spikes./repmat(min(norm_spikes),20,1);

clear pos

x=cov(norm_spikes');
[V,~]=eig(x);
pc_vectors=V(:,end-2:end);
pc1_local=norm_spikes'*pc_vectors(:,end);
pc2_local=norm_spikes'*pc_vectors(:,end-1);
figure
% subplot(2,2,1)
% plot(pc1_local,pc2_local,'.','markersize',1)
% 
isi=round(diff(1000/str2num(timeUnit)*isi_all(rgc{checked}))); % in ms
% 
% for i=0:50
%     a(i+1)=sum(isi==i);
% end
% subplot(2,2,2)
% plot(0:50,a)

% subplot(2,2,3)
a=isi<3&isi>0.1;
b=unique([find(a),find(a)+1]);
plot(pc1_local,tvector,'.','markersize',1)
hold on
plot(pc1_local(b),tvector(b),'.r','markersize',1)

% subplot(2,2,4)
% hist(isi(a),50)
% isi(isi<3)
% plot(pc1(rgc{checked}(a)),1:sum(a),'.','markersize',1)
