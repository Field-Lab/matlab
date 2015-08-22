
% make a list of stimuli.lisp file in data and archive
alldirs=dir('/Volumes/Archive/20*-*-*-*');
cnt=1;
toanalyze=cell(1600,1);
short_name=cell(1600,1);
for i=1:length(alldirs)
    paths=['/Volumes/Archive/',alldirs(i).name,'/Visual/stimuli.lisp'];
    if exist(paths, 'file')
        toanalyze{cnt}=paths;
        cnt=cnt+1;
    end    
end

alldirs=dir('/Volumes/Data/20*-*-*-*');
for i=1:length(alldirs)
    paths=['/Volumes/Data/',alldirs(i).name,'/Visual/stimuli.lisp'];
    if exist(paths, 'file')
        toanalyze{cnt}=paths; 
        cnt=cnt+1;
    end
end

toanalyze(cnt:end)=[];

my_list = cell(900,1);
cnt = 1;
for i=1:length(toanalyze)
    
    fid=fopen(toanalyze{i},'r');
    a=textscan(fid,'%s','Delimiter','\n');
    fclose(fid);
    a = a{1};
    for j=1:length(a)
        if ~isempty(regexpi(a{j}, 'movie_name', 'once'))
            my_list{cnt} = toanalyze{i};
            cnt = cnt+1;
            break;
        end
    end
end