
alldirs=dir('/Volumes/Data/2015-*-*-*');
cnt=1;
for i=1:length(alldirs)
    paths=['/Volumes/Data/',alldirs(i).name,'/Visual/stimuli.lisp'];
    if exist(['/Volumes/Data/',alldirs(i).name,'/Electrical'], 'dir')
        toanalyze{cnt} = [paths, ' Electrical'];
    elseif ~exist(paths, 'file')
        toanalyze{cnt} = [paths,' NO stimuli.lisp'];
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



% fid1=fopen('/Volumes/Analysis/shared/stimuli_db.txt','w');

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

fid1=fopen('/Users/alexth/Desktop/stimuli_db.txt','w');

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
                sizesum=0;
                for kk = 1:length(dataFolder)
                    if dataFolder(kk).isdir
                        subfolders=dir([path2stim,dataFolder.name,'/d*']);                        
                        for j=1:length(subfolders)
                            sizesum=sizesum+subfolders(j).bytes;
                        end
                    else
                        sizesum=dataFolder.bytes;
                    end
                end
            end
            data_sizes(i)=sizesum/1073741824; % in Gb
            
            param_content{overall_cnt,4}=dataFolder.name;
            param_content{overall_cnt,5}=data_sizes(i);
            
            fprintf(fid1, '%s\t%.3f Gb\r\n',...
                dataname(i).name,data_sizes(i));
            
            overall_cnt=overall_cnt+1;
            
        end
        fprintf(fid1, '\r\n\r\n\r\n\r\n');
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
%             elseif ~isempty(regexp(a{i}, '61', 'once')) && isempty(regexp(a{i}, '610', 'once')) &&...
%                     isempty(regexp(a{i}, '; 61 stimulus', 'once')) && isempty(regexp(a{i}, 'data061', 'once'))
%                 err_str=[tmp{1}, '  array 61' ];
%                 %             disp(err_str)
%                 fprintf(fid1, '%s\r\n', 'array 61');
%                 fprintf(fid1, '\r\n\r\n\r\n\r\n');
%                 array61=array61+1;
%                 flagno61=0;
%                 break
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
            dataname=cell(1000,1);datastr=[];
            expr_h='data\d*';
            expr_h1='data \d*';
            expr_miss='_data\d*';
            expr_miss2=':data-*';
            expr_path = 'path';
            cnt=1;
            for i=mystr%comments
                b=regexp(a{i},expr_h1,'match');
                if isempty(b)
                    b=regexp(a{i},expr_h,'match');                    
                end
                c=regexp(a{i},expr_miss,'match');
                d=regexp(a{i},expr_miss2,'match');
                f = regexp(a{i},expr_path,'match');
                if ~isempty(b) && isempty(c) && isempty(d) && isempty(f)                    
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
                        
                        
                    else  % check for grey  or skipped in all lines
                        for j=datastr(i):tstop
                            if ~isempty(regexpi(a{j},'grey|gray|spont|mean'));
                                stim_type='grey';
                                break
                            elseif ~isempty(regexpi(a{j},'skipped'));
                                stim_type='skipped';
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
param_content(overall_cnt:end,:) = [];
comb=[field_names; param_content];


save('/Users/alexth/Desktop/raw_matlab_table_2015.mat','comb')

fid1=fopen('/Users/alexth/Desktop/stimuli_db.txt','w');

for i=1:size(comb,1)
    for j=1:25
        tmp =comb{i,j};
        if iscell(tmp)
            tmp = tmp{1};
        end
        if isnumeric(tmp)
            if isinteger(tmp)
                fprintf(fid1, '%d\t', tmp);
            else
                fprintf(fid1, '%.3f\t', tmp);
            end
        elseif ischar(tmp)
            fprintf(fid1, '%s\t', tmp);
        end
        if j==5
            fprintf(fid1, '\t\t\t\t\t\t\t', tmp);
        end
    end
    fprintf(fid1, '\r\n')
end
fclose(fid1)


