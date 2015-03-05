date='20121004'
codeWord='seq'
path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];
load([path2save,'LF_',codeWord])
% find cells with stable filters at each ND
figure
for cnt=1:size(LinearFilter,4)
    subplot(5,6,cnt)
    hold on
    a=reshape(LinearFilter(:,1,1:24:floor(size(LinearFilter,3)/24)*24,cnt),500,floor(size(LinearFilter,3)/24));
    plot(a)
end 
