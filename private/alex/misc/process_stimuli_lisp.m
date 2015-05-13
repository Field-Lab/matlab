
% make a list of stimuli.lisp file in data and archive
alldirs=dir('/Volumes/Archive/20*-*-*-*');
cnt=1;
toanalyze=cell(1600,1);
short_name=cell(1600,1);
for i=1:length(alldirs)
    paths=['/Volumes/Archive/',alldirs(i).name,'/Visual/stimuli.lisp'];
    if ~exist(paths, 'file')
        toanalyze{cnt} = [paths,' NO stimuli.lisp'];
    elseif exist(['/Volumes/Archive/',alldirs(i).name,'/Electrical'], 'dir')
        toanalyze{cnt} = [paths, ' Electrical'];
    else
        toanalyze{cnt}=paths;
    end
    tmp=regexpi(paths,'\d{4}-\d{2}-\d{2}-\d*','match');
    tmp{1}(regexp(tmp{1},'-'))=[];
    if length(tmp{1})==9
        tmp{1}=[tmp{1}(1:8), '0', tmp{1}(9)];
    end
    short_name{cnt}=tmp{1};
    cnt=cnt+1;
end

alldirs=dir('/Volumes/Data/20*-*-*-*');
for i=1:length(alldirs)
    paths=['/Volumes/Data/',alldirs(i).name,'/Visual/stimuli.lisp'];
    if ~exist(paths, 'file')
        toanalyze{cnt} = [paths,' NO stimuli.lisp'];
    elseif exist(['/Volumes/Data/',alldirs(i).name,'/Electrical'], 'dir')
        toanalyze{cnt} = [paths, ' Electrical'];
    else
        toanalyze{cnt}=paths;        
    end
    tmp=regexpi(paths,'\d{4}-\d{2}-\d{2}-\d*','match');
    tmp{1}(regexp(tmp{1},'-'))=[];
    if length(tmp{1})==9
        tmp{1}=[tmp{1}(1:8), '0', tmp{1}(9)];
    end
    short_name{cnt}=tmp{1};
    cnt=cnt+1;
end

toanalyze(cnt:end)=[];
short_name(cnt:end)=[];

% sort in chronological order
[~,ic]=sort(short_name);
toanalyze=toanalyze(ic);
clear short_name alldirs paths ic tmp



fid1=fopen('/Volumes/Analysis/shared/stimuli_db.txt','w');

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

% exp='(?<=rgb \(mul |rgb #\()(\d\.\d*\s\d\.\d*\s\d\.\d*|\d\.\d*)';
% exp_m='(?<=back-rgb #\()(\d\.\d*\s\d\.\d*\s\d\.\d*|\d\.\d*)';
% exp_t='(?<=[^back-]rgb \(mul |[^back-]rgb #\()(\d\.\d*\s\d\.\d*\s\d\.\d*|\d\.\d*)';
% str='(let ((rgb (mul 0.48 #(1 1 1))))';
% regexpi(str,exp_t,'match')


field_names={'Piece', 'Location', 'exp type','datarun', 'size', 'type',...
    'BW/RGB', 'interval', 'rgb', 'seed', 'back rgb'...
    'stixel width', 'stixel height','field width', 'field height',...
    'x-start','x-end','y-start','y-end','frames','loop','movies',...
    'sparse','jitter', 'double'};
param_content=cell(3000,length(field_names));
shifts=[7, 8, 10, 12, 13, 14, 15, 9, 11, 16, 17, 18, 19, 20, 21, 22, 22, 23, 24];

parsed_cnt=0;
rep_cnt=0;
nodataxxx=0;
array61=0;
wrong_data=0;
overall_cnt=1;
for myDay=length(toanalyze):-1:1
    
    path2data=toanalyze{myDay};
    
    tmp=regexpi(path2data,'\d{4}-\d{2}-\d{2}-\d*','match');
    fprintf(fid1, '%s\r\n',tmp{1});
    
    locs=regexp(path2data,'/');
    locs=path2data(locs(2)+1:locs(3)-1);
    
    param_content(overall_cnt,1)=tmp;
    param_content{overall_cnt,2}=locs;
    
    
    if ~isempty(regexpi(path2data,'NO stimuli.lisp','match'))        
        fprintf(fid1, 'no stimuli.lisp\r\n');
        fprintf(fid1, '\r\n\r\n\r\n\r\n');
        
        param_content{overall_cnt,3}='NO stimuli.lisp';
        overall_cnt=overall_cnt+1;
        
    elseif ~isempty(regexpi(path2data,'Electrical','match'))
        
        param_content{overall_cnt,3}='Electrical';
        
        path2stim=toanalyze{myDay}(1:end-30);
        dataname=dir([path2stim,'data*']);
        data_sizes=zeros(1,length(dataname));
        
        % make template for size
        for i=1:length(dataname)
            dataFolder=dir([path2stim,dataname(i).name,'*']);
           
            if isempty(dataFolder)
                sizesum = 0;
            else
                %look inside folder
                if dataFolder.isdir
                    subfolders=dir([path2stim,dataFolder.name,'/d*']);
                    sizesum=0;
                    for j=1:length(subfolders)
                        sizesum=sizesum+subfolders(j).bytes;
                    end
                else
                    sizesum=dataFolder.bytes;
                end
            end
            data_sizes(i)=sizesum/1073741824; % in Gb
            
            param_content{overall_cnt,4}=dataFolder.name;
            param_content{overall_cnt,5}=data_sizes(i);
            
            overall_cnt=overall_cnt+1;
            
        end
                
    else
                
        param_content{overall_cnt,3}='Visual';
        
        fid=fopen(toanalyze{myDay},'r');
        a=textscan(fid,'%s','Delimiter','\n');
        fclose(fid);
        
        a=a{1, 1};
        
        %remove empty strings
        flagno61=1;
        emptystr=[];
        for i=1:length(a)
            if isempty(a{i})
                emptystr=[emptystr,i];
            elseif ~isempty(regexp(a{i}, '61', 'once')) && isempty(regexp(a{i}, '610', 'once')) && isempty(regexp(a{i}, '; 61 stimulus', 'once'))
                err_str=[tmp{1}, '  array 61' ];
                %             disp(err_str)
                fprintf(fid1, '%s\r\n', 'array 61');
                fprintf(fid1, '\r\n\r\n\r\n\r\n');
                array61=array61+1;
                flagno61=0;
                break
            end
        end
        
        if flagno61
            
            a(emptystr)=[];
            mystr=1:length(a);
            
            
            %find commented strings
            comments=[];
            for i=1:length(a)
                if strcmp(a{i}(1),';')
                    comments=[comments,i];
                end
            end
            
            % find ;data0xx strings (if absent, might crash things!)
            expr_h='data\d*';dataname=cell(50,1);datastr=[];
            expr_miss='_data\d*';
            expr_miss2=':data-*';
            cnt=1;
            for i=mystr%comments
                b=regexp(a{i},expr_h,'match');
                c=regexp(a{i},expr_miss,'match');
                d=regexp(a{i},expr_miss2,'match');
                if ~isempty(b) && isempty(c) && isempty(d)
                    dataname{cnt}=b{:};
                    datastr=[datastr; i];
                    cnt=cnt+1;
                end
            end
            dataname(cnt:end)=[];
            ndata=cnt-1; % number of dataruns
            
            if ndata==0
                err_str=[tmp{1}, '  has no dataXXX!' ];
                disp(err_str)
                fprintf(fid1, '%s\r\n', 'has no dataXXX');
                fprintf(fid1, '\r\n\r\n\r\n\r\n');
                nodataxxx=nodataxxx+1;
                param_content{overall_cnt,3}='no data found in lisp';
                overall_cnt=overall_cnt+1;
            elseif length(unique(dataname))<length(dataname)
                err_str=[tmp{1}, '   has repeated dataXXX!' ];
                disp(err_str)
                fprintf(fid1, '%s\r\n', 'has repeated dataXXX');
                fprintf(fid1, '\r\n\r\n\r\n\r\n');
                rep_cnt=rep_cnt+1;
                param_content{overall_cnt,3}='repeated data ID in lisp';
                overall_cnt=overall_cnt+1;                
            elseif str2num(dataname{end}(end-1:end))~=(length(dataname)-1)
                err_str=[tmp{1}, '   has messed up dataXXX order!' ];
                disp(err_str)
                fprintf(fid1, '%s\r\n', 'wrong numbering of dataXXX');
                fprintf(fid1, '\r\n\r\n\r\n\r\n');
                wrong_data=wrong_data+1;
                param_content{overall_cnt,3}='data IDs wrong in lisp';
                overall_cnt=overall_cnt+1;
            else % parse data
                
                parsed_cnt=parsed_cnt+1;
                uncom=mystr;
                uncom(comments)=[];
                
                path2stim=toanalyze{myDay}(1:end-19);
                mystim=cell(ndata,1);cnt=1;
                data_sizes=zeros(ndata,1);
                
                % identify type and look up size
                for i=1:ndata
                    
                    % size
                    dataFolder=dir([path2stim,dataname{i},'*']);
                    if isempty(dataFolder)
                        sizesum = 0;
                    else
                        %look inside folder
                        if dataFolder.isdir
                            subfolders=dir([path2stim,dataFolder.name,'/d*']);
                            sizesum=0;
                            for j=1:length(subfolders)
                                sizesum=sizesum+subfolders(j).bytes;
                            end
                        else
                            sizesum=dataFolder.bytes;
                        end
                    end
                    data_sizes(i)=sizesum/1073741824; % in Gb
                    
                    param_content(overall_cnt,4)=dataname(i);
                    param_content{overall_cnt,5}=data_sizes(i);
                    overall_cnt=overall_cnt+1;
                    
                    % parse type
                    if i==ndata
                        tstop=length(a);
                    else
                        tstop=datastr(i+1)-1;
                    end
                    
                    stim_type='unknown';
                    if tstop-datastr(i)>3
                        
                        my_string_list=datastr(i)+1:tstop;
                        my_string_list=intersect(my_string_list,uncom);
                        
                        for j=my_string_list
                            mystring=a{j};
                            tmp=regexpi(mystring,expr_type,'match');
                            if ~isempty(tmp)
                                stim_type=tmp{1};
                                break
                            end
                        end
                        
                        if strcmpi(stim_type,'unknown') % if unknown check for grey in all lines
                            for j=datastr(i):tstop
                                if ~isempty(regexpi(a{j},'grey|gray|spont|mean'));
                                    stim_type='grey';
                                    break
                                end
                            end
                        end
                        
                        if strcmpi(stim_type,'binary') % check if repeats
                            for j=my_string_list
                                mystring=a{j};
                                if ~isempty(strfind(mystring,'(repeats '))
                                    stim_type='WN repeats';
                                    break
                                end
                            end
                        end
                        
                        if strcmpi(stim_type,'binary') % check if voronoi WN
                            for j=my_string_list
                                mystring=a{j};
                                if ~isempty(strfind(mystring,'(setq *map* '))
                                    stim_type='Voronoi WN';
                                    break
                                end
                            end
                        end
                        
                        
                    else  % check for grey in all lines
                        for j=datastr(i):tstop
                            if ~isempty(regexpi(a{j},'grey|gray|spont|mean'));
                                stim_type='grey';
                                break
                            end
                        end
                    end
                    
                    mystim{cnt}=stim_type;
                    cnt=cnt+1;
                end
                
                
                
                % analyze stimulus
                params=cell(ndata,22);
                for i=1:ndata
                    if strcmpi(mystim{i}, 'binary') || strcmpi(mystim{i}, 'WN repeats') || strcmpi(mystim{i}, 'movie') 
                        if i==ndata
                            tstop=length(a);
                        else
                            tstop=datastr(i+1);
                        end
                        
                        my_string_list=datastr(i)+1:tstop;
                        my_string_list=intersect(my_string_list,uncom);
                        for j=my_string_list
                            mystring=a{j};
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
                
                % put parameters into table                
                overall_cnt=overall_cnt-ndata;
                for i=1:ndata
                    param_content(overall_cnt,6)=mystim(i);
                    for k = 2:20
                        if length(params{i,k})>1 && k~=9 && k~=10 && k~=17 && k~=18
                            param_content{overall_cnt,shifts(k-1)}=strjoin(params{i,k},', ');
                            param_content{overall_cnt,25}='doubles';
                        elseif ~isempty(params{i,k}) 
                            if k==2
                                if regexpi(params{i,2}{1},'NIL')
                                    param_content{overall_cnt,shifts(k-1)}='BW';
                                else
                                    param_content{overall_cnt,shifts(k-1)}='RGB';
                                end
                            else
                                param_content(overall_cnt,shifts(k-1))=params(i,k);
                            end
                        elseif  ~isempty(params{i,k}) && (k==9 || k==10 || k==17 || k==18)
                            param_content{overall_cnt,shifts(k-1)}=strjoin(params{i,k},', ');
                        end  
                        if k==9 && ~isempty(params{i,k}) && ~isempty(regexp(params{i,k}{1},'\d\.\d*\s\d\.\d*\s\d\.\d*', 'once'))
                            tt=regexp(params{i,k}{1},'\d\.\d*','match');
                            if length(unique(tt))==1
                                param_content{overall_cnt,shifts(k-1)}=unique(tt);
                            end
                        end
                    end
                    if ~isempty(params{i,21}) && ~isempty(params{i,22})
                        param_content{overall_cnt,25}='doubles';
                    end
                    overall_cnt=overall_cnt+1;
                end
                
                % process params into WN movie file name
                frames=cell(ndata,1);
                reps=cell(ndata,1);
                sparsed=cell(ndata,1);
                jitter=cell(ndata,1);
                for i=1:ndata
                    % repeats
                    if ~isempty(params{i,16})
                        reps{i} = params{i,16}{1};
                    else
                        reps{i} = '';
                    end
                    % frames
                    if ~isempty(params{i,15})
                        frames{i} = params{i,15}{1};
                    else
                        frames{i} = '';
                    end
                    % sparse
                    if ~isempty(params{i,19})
                        sparsed{i} = ['sparse! ', params{i,19}{1}];
                    else
                        sparsed{i} = '';
                    end
                    % jitter
                    if ~isempty(params{i,20})
                        jitter{i} = params{i,20}{1};
                    else
                        jitter{i} = '';
                    end
                    if strcmpi(mystim{i}, 'binary') || strcmpi(mystim{i}, 'WN repeats')
                        if ~isempty(find(cellfun(@length, params(i,:))>1, 1))
                            % some double entries?
                            mystim{i,2}='doubles';
                        else
                            % independent?
                            if regexpi(params{i,2}{1},'NIL')
                                myst='BW-';
                            else
                                myst='RGB-';
                            end
                            
                            myst=[myst params{i,5}{1}, '-']; % add stixel size
                            
                            if ~isempty(params{i,3})
                                myst=[myst params{i,3}{1}, '-']; % add refresh rate
                            else
                                myst=[myst 'unk refresh-'];
                                disp([toanalyze{myDay}, ' has no refresh'])
                            end
                            
                            myst=[myst params{i,9}{1}, '-']; % add contrast
                            
                            myst=[myst params{i,4}{1}, '-']; % add seed
                            
                            myst=[myst params{i,7}{1}, 'x']; % add field width
                            
                            myst=[myst params{i,8}{1}]; % add field height
                            
                            mystim{i,2}=myst;
                        end
                    elseif strcmpi(mystim{i}, 'movie')
                        if ~isempty(params{i,17})
                            mystim{i,2}=strjoin(params{i,17},', ');
                        elseif ~isempty(params{i,18})
                            mystim{i,2}=strjoin(params{i,18},', ');
                        else
                            mystim{i,2} = 'movie not found';
                            disp([toanalyze{myDay}, ' has no movie'])
                        end
                    else
                        mystim{i,2}='';
                    end
                end
                
                for i=1:ndata
                    fprintf(fid1, '%s\t%.3f Gb\t%s\t%s\t%s\t%s\t%s\t%s\r\n',...
                        dataname{i},data_sizes(i),frames{i},reps{i}, ...
                        sparsed{i},jitter{i},mystim{i,1}, mystim{i,2});
                end
                
                fprintf(fid1, '\r\n\r\n\r\n\r\n');
            end
        else
            param_content{overall_cnt,3}='Array 61';
            overall_cnt=overall_cnt+1;
        end
    end
end

fclose(fid1)

comb=[field_names; param_content];

save('/Volumes/Analysis/shared/raw_matlab_table.mat','comb')

%%

load('/Volumes/Analysis/shared/raw_matlab_table.mat')

fid=fopen('/Volumes/Analysis/shared/tmp.txt','a')
for i=1:size(comb,1)
    i
    for j=1:size(comb,2)
        x=comb{i,j};
        
        if iscell(x) && length(x)==1
            x=x{1};
        end
        if iscell(x)
            tt=[];
            for k=1:length(x)
                if isnumeric(x{k})
                    x{k}=num2str(x{k});
                end
                tt=[tt, ', ', x{k}];
            end
            x=tt;
        end
        if isnumeric(x)
            x=num2str(x);
        end
        if ~isempty(x)
            if x(1)==','
                x(1)=[];
            end
            if x(1)==' '
                x(1)=[];
            end
        end
        fprintf(fid,'%s\t',x);
    end
    fprintf(fid,'\r\n');
end
fclose(fid)



%% find folders in /Analysis with stim* subfolder/file

a=dir('/Volumes/Analysis/20*');
stimpresent=cell(length(a),2);
cnt=1;
cnt1=0;

for i=1:length(a)
    if a(i).name(1)~='.' && isdir(['/Volumes/Analysis/',a(i).name])
        b=dir(['/Volumes/Analysis/',a(i).name,'/*stim*']);
        if isempty(b)
            b=dir(['/Volumes/Analysis/',a(i).name,'/*Stim*']);
        end
        if ~isempty(b) && b(1).isdir
            b=dir(['/Volumes/Analysis/',a(i).name,'/',b.name]);
        end
        emp=[];
        for j=1:length(b)
            if b(j).name(1)=='.'
                emp=[emp j];
            else
                tmp=find(b(j).name=='.',1,'last');
                if ~isempty(tmp)
                    b(j).name=b(j).name(1:tmp-1);
                end
            end
        end
        b(emp)=[];
        
        if ~isempty(b) %&& ~strcmp(b(1).name(end-3:end),'.mat')
            vis=dir(['/Volumes/Data/',a(i).name,'/Visual']);
            emp=[];
            for j=1:length(vis)
                if vis(j).name(1)=='.'
                    emp=[emp j];
                end
            end
            vis(emp)=[];
            
            if ~isempty(vis)
                flag=0;
                for j=1:length(b)
                    for k=1:length(vis)
                        tmp=find(vis(k).name=='.',1,'last');
                        if ~isempty(tmp)
                            tmp=vis(k).name(1:tmp-1);
                        else
                            tmp=vis(k).name;
                        end
                        if strcmp(tmp, b(j).name)
                            flag=1;
                            break
                        end
                    end
                    if ~flag
                        stimpresent{cnt,2} = 0;
                        cnt1=cnt1+1;
                        
                    end
                end
            end
            stimpresent{cnt,1}=a(i).name;
            cnt=cnt+1;
        end
    end
end
    
stimpresent(cnt:end,:)=[];

for i=1:length(stimpresent)
    if stimpresent{i,2}==0
        disp(stimpresent{i,1});
        tmp=dir(['/Volumes/Analysis/',stimpresent{i,1}]);
        for j=1:length(tmp)
            if tmp(j).name(1)~='.'
                disp(tmp(j).name)
            end
        end
    end
end

%% compare stimuli folders in /Analysis with /Visual in data and archive

a=dir('/Volumes/Analysis/20*');
stimpresent=cell(length(a),2);
cnt=1;
cnt1=0;

for i=1:length(a)
    if a(i).name(1)~='.' && isdir(['/Volumes/Analysis/',a(i).name])
        b=dir(['/Volumes/Analysis/',a(i).name,'/*stim*']);
        if isempty(b)
            b=dir(['/Volumes/Analysis/',a(i).name,'/*Stim*']);
        end
        % get rid of . files
        if ~isempty(b)
            emp=[];
            for j=1:length(b)
                if b(j).name(1)=='.'
                    emp=[emp j];
                else
                    tmp=find(b(j).name=='.',1,'last');
                    if ~isempty(tmp)
                        b(j).name=b(j).name(1:tmp-1);
                    end
                end
            end
            b(emp)=[];
        end
            
        if ~isempty(b) && b(1).isdir % if b is a directory, look inside
            b=dir(['/Volumes/Analysis/',a(i).name,'/',b.name]);
            emp=[];
            for j=1:length(b)
                if b(j).name(1)=='.'
                    emp=[emp j];
                else
                    tmp=find(b(j).name=='.',1,'last');
                    if ~isempty(tmp)
                        b(j).name=b(j).name(1:tmp-1);
                    end
                end
            end
            b(emp)=[];
        end

        
        if ~isempty(b) %&& ~strcmp(b(1).name(end-3:end),'.mat')
            vis=dir(['/Volumes/Archive/',a(i).name,'/Visual']);
            if isempty(vis)
                vis=dir(['/Volumes/Archive/',a(i).name]);
                if ~isempty(vis) % otherwise, do nothing: probably on other volume
                    stimpresent{cnt,1}=a(i).name;
                    stimpresent{cnt,2}='no Visual';
                    stimpresent{cnt,3}=b;
                    cnt=cnt+1;
                end
            else
                emp=[];
                for j=1:length(vis)
                    if vis(j).name(1)=='.'
                        emp=[emp j];
                    else
                        tmp=find(vis(j).name=='.',1,'last');
                        if ~isempty(tmp)
                            vis(j).name=vis(j).name(1:tmp-1);
                        end
                    end
                end
                vis(emp)=[];
                
                if ~isempty(vis)
                    flag=0;
                    for j=1:length(b)
                        for k=1:length(vis)
                            if strcmp(vis(k).name, b(j).name) && b(j).bytes==vis(k).bytes && b(j).datenum==vis(k).datenum
                                flag=1;
                                break
                            end
                        end
                        if ~flag
                            stimpresent{cnt,2} = 0;
                            cnt1=cnt1+1;
                        end
                    end
                end
                stimpresent{cnt,1}=a(i).name;
                cnt=cnt+1;
            end
        end
    end
end
    
stimpresent(cnt:end,:)=[];

for i=1:length(stimpresent)
    if stimpresent{i,2}==0
        disp(stimpresent{i,1});
%         tmp=dir(['/Volumes/Analysis/',stimpresent{i,1}]);
%         for j=1:length(tmp)
%             if tmp(j).name(1)~='.'
%                 disp(tmp(j).name)
%             end
%         end
    end
end

cnt=1;
no_visual=cell(2,2);
for i=1:length(stimpresent)
    if strcmp(stimpresent{i,2},'no Visual')
        no_visual{cnt,1}=stimpresent{i,1};
        
        no_visual{cnt,2}=stimpresent{i,3};
        cnt=cnt+1;
%         tmp=dir(['/Volumes/Analysis/',stimpresent{i,1}]);
%         for j=1:length(tmp)
%             if tmp(j).name(1)~='.'
%                 disp(tmp(j).name)
%             end
%         end
    end
end


a=dir('/Volumes/Archive/20*');
no_visual_archive=cell(2,2);
cnt=1;
for i=1:length(a)
    if isempty(dir(['/Volumes/Archive/',a(i).name,'/Visual']))
        no_visual_archive{cnt,1}=a(i).name;
        b=dir(['/Volumes/Analysis/',a(i).name]);
        emp=[];
        for j=1:length(b)
            if b(j).name(1)=='.'
                emp=[emp j];
            end
        end
        b(emp)=[];
        no_visual_archive{cnt,2}=b;
        cnt=cnt+1;
    end
end

a=dir('/Volumes/Data/20*');
no_visual_data=cell(2,2);
cnt=1;
for i=1:length(a)
    if isempty(dir(['/Volumes/Data/',a(i).name,'/Visual']))
        
        b=dir(['/Volumes/Analysis/',a(i).name]);
        emp=[];
        for j=1:length(b)
            if b(j).name(1)=='.' %|| (length(b(j).name)>3 && (strcmp(b(j).name(1:4), 'data') ...
                    %|| strcmp(b(j).name(1:4), 'stre') ||  strcmp(b(j).name(end-3:end), 'lias')))
                emp=[emp j];
            end
        end
        b(emp)=[];
        if ~isempty(b)            
            no_visual_data{cnt,1}=a(i).name;
            no_visual_data{cnt,2}=b;
            cnt=cnt+1;
        end
    end
end



a=dir('/Volumes/Analysis/20*');
anal_stim=cell(2,2);
cnt=1;
for i=1:length(a)
    b=dir(['/Volumes/Analysis/',a(i).name]);       
    emp=[];
    for j=1:length(b)
        tmp=regexpi(b(j).name, 'stim');
        if ~isempty(tmp) && b(j).name(1)~='.'
            anal_stim{cnt,2}=b(j).name;
            anal_stim{cnt,1}=a(i).name;
            cnt=cnt+1;
        end
    end
end



a=dir('/Volumes/Analysis/20*');
anal_stim=cell(2,3);
cnt=1;
for i=1:length(a)
    b=dir(['/Volumes/Analysis/',a(i).name]);       
    for j=1:length(b)
        tmp=regexpi(b(j).name, 'stim');
        if ~isempty(tmp) && b(j).name(1)~='.'
            anal_stim{cnt,2}=b(j).name;
            anal_stim{cnt,1}=a(i).name;            
            anal_stim{cnt,3}=b(j).isdir;
            cnt=cnt+1;
        end
    end
end

% not a directory (20)
cnt=1;
for i=1:length(anal_stim)
    if anal_stim{i,3}==0
        disp(anal_stim{i,1}) 
        cnt=cnt+1;
    end
end

% directory (2)
for i=1:length(anal_stim)
    if anal_stim{i,3}
        b=dir(['/Volumes/Analysis/',anal_stim{i,1}, '/', anal_stim{i,2}]);
        emp=[];
        for j=1:length(b)
            if b(j).name(1)=='.'
                emp=[emp j];
            end
        end
        b(emp)=[];
        
        if ~isempty(b)
            c=dir(['/Volumes/Archive/',anal_stim{i,1}, '/Visual']);
            if isempty(c)
                c=dir(['/Volumes/Data/',anal_stim{i,1}, '/Visual']);
            end
            if isempty(c)
                disp(anal_stim{i,1})
            end
        end

    end
end


% check and remove stim from analysis 
for i=1:length(anal_stim)
    if anal_stim{i,3} && ~strcmp(anal_stim{i,1}, '2002-10-02-0') && ~strcmp(anal_stim{i,1}, '2009-02-10-0')
        b=dir(['/Volumes/Analysis/',anal_stim{i,1}, '/', anal_stim{i,2}]);
        emp=[];
        for j=1:length(b)
            if b(j).name(1)=='.'
                emp=[emp j];
            end
        end
        b(emp)=[];
        
        if ~isempty(b)            
            paths=['/Volumes/Archive/',anal_stim{i,1}, '/Visual'];
            c=dir(paths);
            if isempty(c)
                paths=['/Volumes/Data/',anal_stim{i,1}, '/Visual'];
                c=dir(paths);
            end
            if isempty(c)
                disp(anal_stim{i,1})
            else
                emp=[];
                for j=1:length(c)
                    if c(j).name(1)=='.'
                        emp=[emp j];
                    end
                end
                c(emp)=[];
            end
            cnt=0;
            for j=1:length(b)
                analysisName=['/Volumes/Analysis/',anal_stim{i,1}, '/', anal_stim{i,2},'/',b(j).name];
                for k=1:length(c)
                    visualName=[paths, '/', c(k).name];
                    t=evalc(['system(''diff ', analysisName, ' ', visualName, ''');']);
                    if isempty(t)
                        cnt=cnt+1;
                    end
                end
            end
            if cnt==length(b) % safe to remove
                [status]=rmdir(['/Volumes/Analysis/',anal_stim{i,1}, '/', anal_stim{i,2}], 's');
            end
            if status==0
                anal_stim{i,1}
            end
        end

    end
end

% get names on Data
a=dir('/Volumes/Data/20*');
names_collection=cell(2,1);
cnt=1;
for i=1:length(a)    
    b=dir(['/Volumes/Data/', a(i).name, '/Visual']);    
    
    emp=[];
    for j=1:length(b)
        if b(j).name(1)=='.' || b(j).isdir || ~isempty(regexpi(b(j).name,'data\w*'))...
                || b(j).name(end)=='w' || (length(b(j).name)>2 && strcmp(b(j).name(end-1:end),'.m'))...
                || b(j).name(end)=='t' || (length(b(j).name)>9 && strcmp(b(j).name(end-7:end),'rawMovie'))...
                || (length(b(j).name)>1 && ~isempty(regexpi(b(j).name(1:2),'s\d')))...
                || ~isempty(regexpi(b(j).name,'map\w*')) || ~isempty(regexpi(b(j).name,'notes\w*'))...
                || ~isempty(regexpi(b(j).name,'experiment\w*'))
            
            emp=[emp j];
        end
    end
    b(emp)=[];
   
    
    if ~isempty(b)        
        for j=1:length(b)
            flag=1;
            for k=1:length(names_collection)
                if strcmp(b(j).name, names_collection{k})
                    flag=0;
                    break
                end
            end
            if flag
                names_collection{cnt}=b(j).name;
                cnt=cnt+1;
            end
        end
    end
end

% rename on data
a=dir('/Volumes/Data/20*');
for i=1:length(a) 
    paths=['/Volumes/Data/', a(i).name, '/Visual/'];
    b=dir(paths);    
    
    if ~isempty(b)
        for j=1:length(b)
            if strcmp(b(j).name, 'stimuli') || strcmp(b(j).name, 'Stimuli.lisp')...
                    || strcmp(b(j).name, 'stim.lisp') || strcmp(b(j).name, 'stimulus')
                
                oldName=[paths b(j).name];
                newName=[paths, 'stimuli.lisp'];
                if ~exist(newName)
                    movefile(oldName, newName);
                    disp([a(i).name, '  ', b(j).name, ' ----> ', 'stimuli.lisp'])
                else
                    disp([a(i).name, '  stimuli.lisp already exists, along with ', b(j).name])
                end
            end
        end
    end
end



% get names on Archive
a=dir('/Volumes/Archive/20*');
names_collection=cell(2,1);
cnt=1;
for i=1:length(a)    
    b=dir(['/Volumes/Archive/', a(i).name, '/Visual']);    
    
    emp=[];
    for j=1:length(b)
        if b(j).name(1)=='.' || b(j).isdir || ~isempty(regexpi(b(j).name,'data\w*'))...
                || b(j).name(end)=='w' || (length(b(j).name)>2 && strcmp(b(j).name(end-1:end),'.m'))...
                || b(j).name(end)=='t' || (length(b(j).name)>9 && strcmp(b(j).name(end-7:end),'rawMovie'))...
                || (length(b(j).name)>1 && ~isempty(regexpi(b(j).name(1:2),'s\d')))...
                || ~isempty(regexpi(b(j).name,'map\w*')) || ~isempty(regexpi(b(j).name,'notes\w*'))...
                || ~isempty(regexpi(b(j).name,'experiment\w*'))...
                || (length(b(j).name)>4 && ~isempty(regexpi(b(j).name(end-3:end),'.xls')))...
                || ~isempty(regexpi(b(j).name,'mb\w*'))
            
            emp=[emp j];
        end
    end
    b(emp)=[];
   
    
    if ~isempty(b)        
        for j=1:length(b)
            flag=1;
            for k=1:length(names_collection)
                if strcmp(b(j).name, names_collection{k})
                    flag=0;
                    break
                end
            end
            if flag
                names_collection{cnt}=b(j).name;
                cnt=cnt+1;
            end
        end
    end
end

% rename on Archive
a=dir('/Volumes/Archive/20*');
for i=1:length(a) 
    paths=['/Volumes/Archive/', a(i).name, '/Visual/'];
    b=dir(paths);    
    
    if ~isempty(b)
        for j=1:length(b)
            if strcmp(b(j).name, 'stimuli') || strcmp(b(j).name, 'Stimuli.lisp')...
                    || strcmp(b(j).name, 'stim.lisp') || strcmp(b(j).name, 'stimulus')...
                    || strcmp(b(j).name, 'stim.lisp') || strcmp(b(j).name, 'Stimulus.lisp')...
                    || strcmp(b(j).name, 'stiimuli.lisp') || strcmp(b(j).name, 'stiimulus.lisp')...
                    || strcmp(b(j).name, 'stim') || strcmp(b(j).name, 'stimulus.lisp')
                
                oldName=[paths b(j).name];
                newName=[paths, 'stimuli.lisp'];
                if ~exist(newName)
                    movefile(oldName, newName);
                    disp([a(i).name, '  ', b(j).name, ' ----> ', 'stimuli.lisp'])
                else
                    disp([a(i).name, '  stimuli.lisp already exists, along with ', b(j).name])
                end
            end
        end
    end
end

% search for specific name on Archive
a=dir('/Volumes/Archive/20*');
for i=1:length(a) 
    paths=['/Volumes/Archive/', a(i).name, '/Visual/'];
    b=dir(paths);    
    
    if ~isempty(b)
        for j=1:length(b)
            if strcmp(b(j).name, 'clouds-temporary.lisp') 
                disp(a(i).name)
                break
            end
        end
    end
end



% check if there is /Visual on Archive/Data (change it) which has no stimuli.lisp
a=dir('/Volumes/Data/20*');
cnt=1;
for i=1:length(a) 
    if exist(['/Volumes/Data/', a(i).name, '/Visual/'], 'dir')
        if ~exist(['/Volumes/Data/', a(i).name, '/Visual/stimuli.lisp'], 'file')
            disp(a(i).name)
            cnt=cnt+1;
        end
    end
end


%%
names_collection=cell(2,1);
cnt=1;
for i=1:length(no_visual)    
    tmp=no_visual{i,2};
    for j=1:length(tmp)
        flag=1;
        for k=1:length(names_collection)
            if strcmp(tmp(j).name, names_collection{k})
                flag=0;
                break
            end
        end
        if flag
             names_collection{cnt}=tmp(j).name;
             cnt=cnt+1;
        end
    end
end

% find files with only typical visual stimuli
my_visual=[];
for i=1:length(no_visual)    
    tmp=no_visual{i,2};
    jj=0;
    for j=1:length(tmp)
        if ~isempty(regexpi(tmp(j).name,'s\d*|stimuli|map'))
            jj=jj+1;            
        end
    end
    if jj==length(tmp)
        my_visual=[my_visual; no_visual{i,1}];
    end
end
% move them

for i=1:size(my_visual,1)
    if isdir(['/Volumes/Archive/',my_visual(i,:),'/Visual'])
        my_visual(i,:)
    else
        tmp=dir(['/Volumes/Analysis/',my_visual(i,:),'/Stimuli/']);
        pths=['/Volumes/Analysis/',my_visual(i,:),'/Stimuli/'];
        if isempty(tmp)
            tmp=dir(['/Volumes/Analysis/',my_visual(i,:),'/stimuli/']);
            pths=['/Volumes/Analysis/',my_visual(i,:),'/stimuli/'];
            if isempty(tmp)
                disp(my_visual(i,:))
                disp('No analysis folder')
            end
        end
                
        if ~isempty(tmp)
            emp=[];
            for j=1:length(tmp)
                if tmp(j).name(1)=='.'
                    emp=[emp j];
                end
            end
            tmp(emp)=[];
            
            if ~isempty(tmp)                
                if ~exist(['/Volumes/Archive/',my_visual(i,:),'/Visual/'],'dir')
                    mkdir(['/Volumes/Archive/',my_visual(i,:),'/Visual/']);
                end
            
                for j=1:length(tmp) 
                    if strcmp(tmp(j).name, 'stim.lisp') || strcmp(tmp(j).name, 'stim') ...
                            || strcmp(tmp(j).name, 'stimulus.lisp') || strcmp(tmp(j).name, 'stimulus')
                        disp(my_visual(i,:))
                        disp([tmp(j).name ' to stimuli.lisp'])
                        newName=['/Volumes/Archive/',my_visual(i,:),'/Visual/stimuli.lisp'];
                    else
                        newName=['/Volumes/Archive/',my_visual(i,:),'/Visual/' tmp(j).name];
                    end
                    oldName=[pths tmp(j).name];       
                    if ~exist(newName)
                        movefile(oldName, newName);
                    else
                        disp( my_visual(i,:))
                        disp('CHECK!!!')
                        
                    end
                end
            end
            tmp=dir(pths);
            emp=[];
            for j=1:length(tmp)
                if tmp(j).name(1)=='.'
                    emp=[emp j];
                end
            end
            tmp(emp)=[];
            if ~isempty(tmp)
                my_visual(i,:)
                tmp
            else
                rmdir(pths)
            end
        end
    end
end



% find files with electrical stims
cnt=1;
my_electric=[];
for i=1:length(no_visual)    
    tmp=no_visual{i,2};
    for j=1:length(tmp)
        if ~isempty(regexpi(tmp(j).name,'pattern\d{3}|movie\d{3}'))
            my_electric=[my_electric; no_visual{i,1}];
            break
        end
    end
end

for i=1:size(my_electric,1)
    if isdir(['/Volumes/Archive/',my_electric(i,:),'/Electrical']) || isdir(['/Volumes/Archive/',my_electric(i,:),'/Visual'])
        my_electric(i,:)
    else
        tmp=dir(['/Volumes/Analysis/',my_electric(i,:),'/stim_files/']);
        if isempty(tmp)
            my_electric(i,:)
        else
            emp=[];
            for j=1:length(tmp)
                if tmp(j).name(1)=='.'
                    emp=[emp j];
                end
            end
            tmp(emp)=[];
            for j=1:length(tmp)
                if ~isempty(regexpi(tmp(j).name, 'pattern\d{3}|movie\d{3}|inputs'))
                    if ~exist(['/Volumes/Archive/',my_electric(i,:),'/Electrical/'],'dir')
                        mkdir(['/Volumes/Archive/',my_electric(i,:),'/Electrical/']);
                    end
                    
                    oldName=['/Volumes/Analysis/',my_electric(i,:),'/stim_files/' tmp(j).name];
                    newName=['/Volumes/Archive/',my_electric(i,:),'/Electrical/'  tmp(j).name];
                    movefile(oldName, newName);
                else
                    if ~exist(['/Volumes/Archive/',my_electric(i,:),'/Visual/'],'dir')
                        mkdir(['/Volumes/Archive/',my_electric(i,:),'/Visual/']);
                    end
                    if ~isempty(regexpi(tmp(j).name, '\w*.lisp'))
                        oldName=['/Volumes/Analysis/',my_electric(i,:),'/stim_files/' tmp(j).name];
                        newName=['/Volumes/Archive/',my_electric(i,:),'/Visual/stimuli.lisp'];
                        movefile(oldName, newName);
                    end
                end
            end
            tmp=dir(['/Volumes/Analysis/',my_electric(i,:),'/stim_files/']);
            emp=[];
            for j=1:length(tmp)
                if tmp(j).name(1)=='.'
                    emp=[emp j];
                end
            end
            tmp(emp)=[];
            if ~isempty(tmp)
                my_electric(i,:)
                tmp
            else
                rmdir(['/Volumes/Analysis/',my_electric(i,:),'/stim_files/'])
            end
        end
    end
end

%% rename fodlers without piece number to -0

a=dir('/Volumes/Data/20*');
cnt=1;
noperm=cell(2,2);
for i=1:length(a)
    if isdir(['/Volumes/Data/',a(i).name])
        if length(a(i).name)==10 % no piece number            
            b=dir(['/Volumes/Data/', a(i).name, '-*']); % check if there are folders for same date with piece numbers
            if isempty(b) % nothing, rename assigning piece 0
                oldName=['/Volumes/Data/' a(i).name];
                newName=[oldName, '-0'];
                try
                    movefile(oldName, newName);
                catch err
                    noperm{cnt,1}=oldName;
                    noperm{cnt,2}=err;
                    cnt=cnt+1;
                end
            else
                a(i).name
            end
        end
    end
end
    

%% find folders in /Data or /Archive which have no /Visual subfolder but have stim* subfolder/file in /Analysis

a=dir('/Volumes/Data/*');
noVis=cell(length(a),1);

cnt=1;
for i=1:length(a)
    if a(i).name(1)~='.' && isdir(['/Volumes/Data/',a(i).name])
        b=dir(['/Volumes/Data/',a(i).name,'/Visual']);
        if isempty(b)
            b=dir(['/Volumes/Analysis/',a(i).name,'/stim*']); 
            if isempty(b)
                b=dir(['/Volumes/Analysis/',a(i).name,'/Stim*']);
            end
            if ~isempty(b)
                noVis{cnt}=a(i).name;               
                cnt=cnt+1;
            end
        end
    end
end

noVis(cnt:end)=[];

%%


alldirs=dir('/Volumes/Data/20*-*-*-*');
cnt=1;cnt1=1;
toanalyze=cell(600,1);nostim=cell(2,3);
for i=1:length(alldirs)
    tmp=dir(['/Volumes/Data/',alldirs(i).name,'/Visual/stim*.lisp']);
    if isempty(tmp)
        tmp=dir(['/Volumes/Data/',alldirs(i).name,'/Visual/',alldirs(i).name(1:10),'*.lisp']);
    end
    if ~isempty(tmp)
        toanalyze{cnt}=['/Volumes/Data/',alldirs(i).name,'/Visual/',tmp.name];
        cnt=cnt+1;
    else    % no stimulus file..?    
        nostim{cnt1,1}=alldirs(i).name;
        % check notebook
        tmp=dir(['/Volumes/Lab/Notebooks/',alldirs(i).name(1:10),'*']);
        if isempty(tmp)
            nostim{cnt1,2}='no notebook';
        end
        % check data files
        tmp=dir(['/Volumes/Data/',alldirs(i).name,'/data*']);
        if isempty(tmp)
            tmp=dir(['/Volumes/Data/',alldirs(i).name,'/*maging*']);
            if isempty(tmp)
                nostim{cnt1,3}='no data folders';
            else
                nostim{cnt1,3}='Imaging only';
            end
        else
            flag=2;
            for j=1:length(tmp)
                tmp2=dir(['/Volumes/Data/',alldirs(i).name,'/',tmp(j).name,'/data*']);
                if ~isempty(tmp2)
                    flag=0;
                    break
                end
            end
            if flag
                
                tmp=dir(['/Volumes/Data/',alldirs(i).name,'/*maging*']);
                if ~isempty(tmp)
                    nostim{cnt1,3}='Imaging only';
                else
                    nostim{cnt1,3}='no data files';
                end
            end
        end
        cnt1=cnt1+1;
    end
end
toanalyze(cnt:end)=[];




% rename files
cnt=1;noperm=cell(3,2);
for i=1:length(toanalyze)

   oldName=toanalyze{i};
   tmp=regexpi(oldName,'/');   
   newName=[oldName(1:tmp(end)),'stimuli.lisp'];
   if ~strcmp(oldName,newName)
       try 
           movefile(oldName, newName);
       catch err  
           noperm{cnt,1}=oldName;
           noperm{cnt,2}=err;
           cnt=cnt+1;
       end
   end
end



alldirs=dir('/Volumes/Archive/20*');
cnt=1;
for i=1:length(alldirs)
    tmp=dir(['/Volumes/Lab/Notebooks/',alldirs(i).name(1:10),'*']);
    if isempty(tmp)
        disp(alldirs(i).name)
        cnt=cnt+1;
    end
end



cnt=1;no_notebook=cell(3,1);
for i=1:length(nostim)
    
    if strcmp(nostim{i,2}, 'no notebook')
        no_notebook{cnt}=nostim{i,1};
        cnt=cnt+1;
    end
end




cnt=1;no_notebook=cell(3,1);
for i=1:length(nostim)
    
    if strcmp(nostim{i,3}, 'no notebook')
        no_notebook{cnt}=nostim{i,1};
        cnt=cnt+1;
    end
end



cnt=1;no_df=cell(3,3);
for i=1:length(nostim)
    
    if strcmp(nostim{i,3}, 'no data folders')
        no_df{cnt,1}=nostim{i,1};
        tmp=dir(['/Volumes/Data/',no_df{cnt},'/']);
        a=ones(size(tmp));
        for j=1:length(tmp)
            if tmp(j).name(1)=='.';
                a(j)=0;
            end
        end
        tmp=tmp(find(a));
        no_df{cnt,3}=tmp(1).name;
        no_df{cnt,2}=length(tmp);
        cnt=cnt+1;
    end
end




cnt=1;no_dfs=cell(3,3);
for i=1:length(nostim)
    
    if strcmp(nostim{i,3}, 'no data files')
        no_dfs{cnt,1}=nostim{i,1};
        tmp=dir(['/Volumes/Data/',no_dfs{cnt},'/']);
        a=ones(size(tmp));
        tt=0;
        for j=1:length(tmp)
            if tmp(j).name(1)=='.';
                a(j)=0;
            else
                tt=tt+tmp(j).bytes;
            end
        end
        tmp=tmp(find(a));
        no_dfs{cnt,3}=tmp(1).name;
        no_dfs{cnt,2}=length(tmp);
        no_dfs{cnt,4}=tt/1024/1000/1000; % in GB
        cnt=cnt+1;
    end
end


