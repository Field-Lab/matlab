function [cones, cone_peaks] = find_cones_auto(sta, threshold)

field_size = size(sta);

if length(field_size)>2 
    sta = sum(sta,3);
    field_size = field_size(1:2);
end

tmp = sta;
marks = threshold * robust_std(tmp(:));

tmp(tmp<marks) = 0;

[k,c] = sort(tmp(:));

[i,j] = ind2sub(field_size, c(find(k>marks,1):end));

tm = find(c>marks*2,1);

ti = sum(c(end-tm:end).* i(end-tm:end))/sum(c(end-tm:end));

tj = sum(c(end-tm:end).* j(end-tm:end))/sum(c(end-tm:end));

a = pdist2([ti,tj], [i,j])';

% figure;
% plot(a)

m = k(end-length(a)+1:end)./(a+1);

% figure
% plot(m)

[~, t] = sort(m);
tt = c(end-length(a)+1:end);

my_indices = tt(t);
% figure
% plot(tmp(my_indices))

window_length = min(10, round(length(t)/2));
moving_window = ones(1, window_length);

a = conv(tmp(my_indices), moving_window, 'same');
k=find(a>max(a(1:round(end/3))),1);

my_indices = my_indices(k:end);
[i,j] = ind2sub(field_size, my_indices);

cones = [j,i];

% hold on
% plot(j,i,'xr')

tmp = zeros(field_size(1), field_size(2));
for k=1:length(i)
    tmp(i(k),j(k)) = sta(i(k),j(k));
end
% figure
% imagesc(tmp)

cone_peaks = find_local_maxima(tmp, 'radius', 2, 'thresh', 0.01, 'return', 'indices');
cone_peaks = [cone_peaks(:,2),cone_peaks(:,1)];


% 
% figure
% colormap gray
% imagesc(sta)
% hold on
% plot(j,i,'xr')
% plot(cone_peaks(:,1),cone_peaks(:,2),'+y')

scale = 5;

for k=1:size(cone_peaks,1)
    tt = sta(cone_peaks(k,2)-1:cone_peaks(k,2)+1,cone_peaks(k,1)-1:cone_peaks(k,1)+1);
%     figure
%     colormap gray
%     imagesc(tt)

    tt = imresize(tt, scale);
    
    [p, m]= find(tt == max(tt(:)),1);
%     [p, m] = [p, m]-8;
%     hold on
%     plot(m/(scale-1), p/(scale-1), 'rx')
    
    cone_peaks(k,1) = cone_peaks(k,1)+(m-8)/scale;
    cone_peaks(k,2) = cone_peaks(k,2)+(p-8)/scale;
%     figure
%     colormap gray
%     imagesc(tt)
end

% plot(cone_peaks(:,1),cone_peaks(:,2),'*g')

tmp = (pdist2(cone_peaks,cone_peaks));
tmp(tmp==0)=100;
p = min(tmp);
t = find(p<2);
m = zeros(1,length(t));
for i=1:length(t)
    m(i) = find(tmp(:,t(i))<2,1);
end
t = unique(sort([t; m])', 'rows');
cone_peaks(t(:,1),:) = [];

