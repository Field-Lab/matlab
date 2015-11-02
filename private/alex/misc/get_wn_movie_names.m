function [my_wn_movie, stixel_size] = get_wn_movie_names(filepath)


fid=fopen(filepath,'r');
my_text=textscan(fid,'%s','Delimiter','\n');
fclose(fid);

my_text=my_text{1, 1};

expr_type='(?<=:type :)(\w*[-]\w*|\w*)';

expr{1}='(?<=:type :)(binary|gaussian)';
expr{2}='(?<=:independent )(NIL|t)';
expr{3}='(?<=:interval )\d*';
expr{4}='(?<=:seed |fixed-seed |seed )\d*';
expr{5}='(?<=:stixel-width )\d*';
expr{6}='(?<=:stixel-height )\d*';
expr{7}='(?<=:field-width )\d*';
expr{8}='(?<=:field-height )\d*';
expr{9}='(?<=[^back-]rgb \(mul |[^back-]rgb #\()(\d\.\d*\s\d\.\d*\s\d\.\d*|\d\.\d*|\d)';
expr{10}='(?<=back-rgb #\()(\d\.\d*\s\d\.\d*\s\d\.\d*|\d\.\d*|\d)';
expr{11}='(?<=:x-start )\d*';
expr{12}='(?<=:x-end )\d*';
expr{13}='(?<=:y-start )\d*';
expr{14}='(?<=:y-end )\d*';
expr{15}='(?<=frames )\d*';
expr{16}='(?<=repeats )\d*';
expr{17}='(?<=path ")\S*.\S*(?=")';
expr{18}='(?<=")\S*.\S*(?=.rawMovie")'; % list of movies
expr{19}='(?<=:probability )(\d*\.\d*|\d*)';
expr{20}='(?<=:jitter )\w*';
expr{21}='(?<=Stimulus-)A';
expr{22}='(?<=Stimulus-)B';


%remove empty strings
emptystr=[];
for i=1:length(my_text)
    if isempty(my_text{i})
        emptystr=[emptystr,i];
    end
end
my_text(emptystr)=[];
mystr=1:length(my_text);

%find dataxxx in commented strings
comments=[];
expr_h='data\d*';dataname=cell(100,1);datastr=[];
expr_miss='_data\d*';
expr_miss2=':data-*';
cnt = 1;
for i=1:length(my_text)
    if strcmp(my_text{i}(1),';')
        comments=[comments,i];
        
        b=regexp(my_text{i},expr_h,'match');
        c=regexp(my_text{i},expr_miss,'match');
        d=regexp(my_text{i},expr_miss2,'match');
        if ~isempty(b) && isempty(c) && isempty(d)
            dataname{cnt}=b{:};
            datastr=[datastr; i];
            cnt=cnt+1;
        end
        
    end
end
dataname(cnt:end)=[];
ndata=cnt-1; % number of dataruns

if ndata==0 || length(unique(dataname))<length(dataname) || str2num(dataname{end}(end-1:end))~=(length(dataname)-1)
    disp('\n\n data numbers incorrect\n\n');
    return
end


uncom=mystr;
uncom(comments)=[];

stim_types = cell(1,ndata);
% identify type and look up size
for i=1:ndata
    
    % parse type
    if i==ndata
        tstop=length(my_text);
    else
        tstop=datastr(i+1)-1;
    end
    
    stim_type='unknown';
    if tstop-datastr(i)>3
        
        my_string_list=datastr(i)+1:tstop;
        my_string_list=intersect(my_string_list,uncom);
        
        for j=my_string_list
            mystring=my_text{j};
            tmp=regexpi(mystring,expr_type,'match');
            if ~isempty(tmp)
                stim_type=tmp{1};
                break
            end
        end
        
        if strcmpi(stim_type,'unknown') % if unknown check for grey in all lines
            for j=datastr(i):tstop
                if ~isempty(regexpi(my_text{j},'grey|gray|spont|mean'));
                    stim_type='grey';
                    break
                end
            end
        end
        
        if strcmpi(stim_type,'binary') % check if repeats
            for j=my_string_list
                mystring=my_text{j};
                if ~isempty(strfind(mystring,'(repeats '))
                    stim_type='WN repeats';
                    break
                end
            end
        end
        
        if strcmpi(stim_type,'binary') % check if voronoi WN
            for j=my_string_list
                mystring=my_text{j};
                if ~isempty(strfind(mystring,'(setq *map* '))
                    stim_type='Voronoi WN';
                    break
                end
            end
        end
        
        
    else  % check for grey in all lines
        for j=datastr(i):tstop
            if ~isempty(regexpi(my_text{j},'grey|gray|spont|mean'));
                stim_type='grey';
                break
            end
        end
    end
    stim_types{i}=stim_type;
    
end


% analyze stimulus
params=cell(ndata,22);
for i=1:ndata
    if strcmpi(stim_types{i}, 'binary') || strcmpi(stim_types{i}, 'WN repeats')...
            || strcmpi(stim_types{i}, 'movie') || strcmpi(stim_types{i}, 'Voronoi WN')
        if i==ndata
            tstop=length(my_text);
        else
            tstop=datastr(i+1);
        end
        
        my_string_list=datastr(i)+1:tstop;
        my_string_list=intersect(my_string_list,uncom);
        for j=my_string_list
            mystring=my_text{j};
            for k=1:22
                b=regexpi(mystring,expr{k},'match');
                if ~isempty(b)
                    if isempty(params{i,k})
                        params{i,k}=b;
                    else
                        params{i,k}=[params{i,k} b];
                    end
                    
                end
            end
        end
        
    end
end

% make white noise movies
my_wn_movie = cell(1,ndata);
stixel_size = zeros(1,ndata);
for i=1:ndata
    if strcmpi(stim_types{i}, 'binary') || strcmpi(stim_types{i}, 'Voronoi WN')% include voronoi in future 
        if strcmp(params{i,2}, 'NIL')
            wnm = 'BW-';
        else
            wnm = 'RGB-';
        end
        
        if ~isempty(params{i,19})
            wnm = [wnm 'sparse'];
        end
        
        if str2num(params{i,7}{1})*str2num(params{i,5}{1})==640 && str2num(params{i,8}{1})*str2num(params{i,5}{1})==320 % don't attach size ending
            wnm = [wnm params{i,5} '-' params{i,3} '-' params{i,9} '-' params{i,4} '.xml']; 
        else        
            wnm = [wnm params{i,5} '-' params{i,3} '-' params{i,9} '-' params{i,4} '-' params{i,7} 'x' params{i,8} '.xml'];
        end
        my_wn_movie{i} = cell2mat(wnm);  
        stixel_size(i) = str2num(params{i,5}{1});
    end
end

