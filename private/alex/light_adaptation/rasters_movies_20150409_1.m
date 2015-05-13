data='011';
datarun = load_data(['/Volumes/Analysis/2015-04-09-1/d11-20-norefit/data',data,'-from-d11_20/data',data,'-from-d11_20']);
datarun = load_params(datarun,'verbose',1);
datarun = set_polarities(datarun);
datarun = load_sta(datarun,'load_sta','all','keep_java_sta',true);
datarun = load_neurons(datarun);
datarun = load_ei(datarun, 'all');

data='011';
datarun1 = load_data(['/Volumes/Analysis/2015-04-09-1/d02-11-norefit/data',data,'-from-d02_11/data',data,'-from-d02_11']);
datarun1 = load_params(datarun1,'verbose',1);
datarun1 = set_polarities(datarun1);
datarun1 = load_sta(datarun1,'load_sta','all','keep_java_sta',true);
datarun1 = load_neurons(datarun1);
datarun1 = load_ei(datarun1, 'all');

[a,b] = map_ei(datarun,datarun1);

figure
for i=1:5
    subplot(2,3,i)
    plot_rf_summaries(datarun, {i}, 'clear', false,  'plot_fits', true, 'fit_color', 'r')
    hold on
    plot_rf_summaries(datarun1, {i}, 'clear', false,  'plot_fits', true, 'fit_color', 'b')
    title(datarun.cell_types{i}.name)
end

sta = zeros(16,16,3,length(datarun.cell_ids));
for i = 1:length(datarun.cell_ids)
    sta(:,:,:,i) = double(datarun.stas.stas{i}(:,:,:,25));
    tt=sta(:,:,:,i);
    sta(:,:,:,i) = sta(:,:,:,i) / max(abs(tt(:)));
end

sta1 = zeros(16,16,3,length(datarun1.cell_ids));
for i = 1:length(datarun1.cell_ids)
    sta1(:,:,:,i) = double(datarun1.stas.stas{i}(:,:,:,25));   
    tt=sta1(:,:,:,i);
    sta1(:,:,:,i) = sta1(:,:,:,i) / max(abs(tt(:)));
end

tmp1=sta(:,:,:,find(datarun.cell_ids==888));
tmp2=sta1(:,:,:,find(datarun1.cell_ids==5896));
tmp3=sta1(:,:,:,find(datarun1.cell_ids==888));

sum((tmp1(:)-tmp2(:)))

figure
imagesc(tmp1(:,:,2))

figure
imagesc(tmp2(:,:,2))

figure
imagesc(tmp3(:,:,2))

sum((tmp1(:)-tmp3(:)))
tt = tmp1(:,:,2,25) - tmp3(:,:,2,25);
sum(abs(tt(:)))


for i = 1:length(datarun.cell_ids)
    tmp = repmat(sta(:,:,:,i), 1,1,1,size(sta1,4)) - sta1;
    tmp = reshape(tmp,numel(tmp)/size(tmp,4),size(tmp,4));
    [sta_diff(i),ind(i)]=min(sum(abs(tmp)));
end

figure
plot(sort(sta_diff))
tmp = sort(sta_diff);
a_sta = ind;
for i=1:length(a)
    if isempty(a{i})
        a_ei(i) = 0;
    else
        a_ei(i) = find(datarun1.cell_ids==a{i})
    end
end

figure
plot(a_ei,a_sta,'*')

tmp=find(a_ei-a_sta~=0)
clear bc
for i=1:length(tmp)
    if a_ei(tmp(i))
        bc(i)=datarun1.cell_ids(a_ei(tmp(i)));
    else
        bc(i) = 0;
    end
end
aa=[datarun.cell_ids(tmp)' datarun1.cell_ids(a_sta(tmp))' bc' sta_diff(tmp)']










ncells=length(datarun.cell_ids);

figure
hold on
g = [];
ss=[];extr=[];
for i=1:ncells
    plot(datarun.vision.timecourses(i).r,'r')    
    plot(datarun.vision.timecourses(i).g,'g')
    plot(datarun.vision.timecourses(i).b,'b')
    g = [g datarun.vision.timecourses(i).g];
    if ~isempty(datarun.vision.timecourses(i).g)
        ss=[ss mean(datarun.vision.sta_fits{i}.sd)];
        extr = [extr datarun.vision.timecourses(i).g(26)];
    end
end

figure
plot(datarun.vision.timecourses(1).g,'g')
hold on
plot(datarun.vision.timecourses(1).b,'b')
plot(datarun.vision.timecourses(1).r,'r')  


g=cov(g);
[V,~]=eig(g);
pc_vectors=V(:,end-2:end);
pc1=g*pc_vectors(:,end);
pc2=g*pc_vectors(:,end-1);
figure
plot(pc1,pc2,'*')

figure
plot(ss,pc1,'*')

figure
plot(ss,extr,'*')

% NDF4
% coarse sta data002
% NSEM data003
% repeats data004
% fine sta data005

% NDF3
% coarse sta data006
% NSEM data021
% repeats data022
% fine sta data020

% NDF2
% coarse sta data007
% NSEM data008
% repeats data009
% fine sta data010

% NDF1
% coarse sta data019
% NSEM data017, data023
% repeats data018
% fine sta data016

% NDF0
% coarse sta data011
% NSEM data012
% repeats data013
% fine sta data014
