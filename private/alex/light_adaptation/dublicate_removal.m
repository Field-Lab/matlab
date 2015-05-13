
% paths

mainpath = '/Volumes/Analysis/2015-03-09-2/';

subpath{1} = 'd02-11-norefit/';
subpath{2} = 'd10-19-norefit/';
subpath{3} = 'd18-28-norefit/';
subpath{4} = 'd27-33-norefit/';
subpath{5} = 'd32-37-norefit/';
subpath{6} = 'd37-40-norefit/';



datarun = load_data([mainpath, subpath{6}, '/data038-from-d37-40/data038-from-d37-40']);
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun,'load_sta',[],'keep_java_sta',true);
datarun = set_polarities(datarun);
datarun = load_neurons(datarun);
datarun = load_ei(datarun,'all');

datarun2 = load_data([mainpath, subpath{2}, '/data011-from-d10-d19/data011-from-d10-d19']);
datarun2 = load_params(datarun2,'verbose',1);
datarun2 = load_sta(datarun2,'load_sta',[],'keep_java_sta',true);
datarun2 = set_polarities(datarun2);
datarun2 = load_neurons(datarun2);
datarun2 = load_ei(datarun2,'all');

mean_sum=zeros(296,296);
i=158;
j=159;
j=156;
for i=1:length(datarun.cell_ids)-1
    for j=i+1:length(datarun.cell_ids)
        a=datarun.ei.eis{i}/min(datarun.ei.eis{i}(:));
        b=datarun.ei.eis{j}/min(datarun.ei.eis{j}(:));
        t=a-b;
        maxx(i,j)=max(t(:));
        minn(i,j)=min(t(:));
    end
end


t=maxx-minn;
imagesc(t)
hist(t(t<0.8&t>=0.7),20)
[a,b]=find(t<0.8&t>0)
tm=[];
for i=1:68
    tm=[tm t(a(i),b(i))];
end

pairs=[a b];

pairs=[(1:68)' a b datarun.cell_ids(a)' datarun.cell_ids(b)' tm'*100 ];

reals2=[2031 2041]

reals=[843 844;2267 2270;2642 2644;2716 2717;3006 3008;3018 3021;3018 3032;...
    3587 3589; 3587 4085; 3589 4085; 4306 4307;4381 4382;4381 4383;4382 4383;...
    4381 4387; 4382 4387;4383 4387;4381 4388;4382 4388;4383 4388; 4387 4388;...
    4388 4478;4937 4938; 4966 4967; 5191 5192;6482 6483; 6513 6663; 6706 6707;...
    6886 6887; 7217 7218;]

maybes=[2266 2272;4111 4115; 5627 5629; 7158 7293]
  
nn=[2 7 11 12 13 14 15 17 18 19 21 22 23 24 25 26 27 28 29 30 31 32 34 35 36 38 39 40 41 42];
tm=[];
for i=1:30
    tm=[tm t(a(nn(i)),b(nn(i)))];
end
figure
hist(tm)
clear amp
for cnt=1:size(pairs,1)
    i=pairs(cnt,2);
    j=pairs(cnt,3);
    tmp=corr(datarun.ei.eis{i}',datarun.ei.eis{j}');
    amp(cnt)=length(tmp(tmp>0.99));
end
plot(tm,amp,'*')
% imagesc(tmp)



amp=zeros(296,296);

for i=1:length(datarun.cell_ids)-1
    i
    for j=i+1:length(datarun.cell_ids)
        tmp=corr(datarun.ei.eis{i}',datarun.ei.eis{j}');
        amp(i,j)=length(tmp(tmp>0.99));
    end
end

hist(amp(:),100)


