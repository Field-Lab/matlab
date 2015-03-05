path2data='/Users/alexth/Desktop/old_stuff/20130703a/20130703/';
a=dir([path2data,'*.mat'])

dest_path1='/Users/alexth/Desktop/old_stuff/20130703a/timing/';
dest_path2='/Users/alexth/Desktop/old_stuff/20130703b/timing/';
if ~exist(dest_path1,'dir')
    mkdir(dest_path1);
end
if ~exist(dest_path2,'dir')
    mkdir(dest_path2);
end

for i=1:length(a)
    k=a(i).name;
    if strcmp(k(end-5:end-4),'_1')
        copyfile([path2data,a(i).name],[dest_path2,k(1:end-6),'.mat'])
    else
        copyfile([path2data,a(i).name],[dest_path1,k])
    end
end