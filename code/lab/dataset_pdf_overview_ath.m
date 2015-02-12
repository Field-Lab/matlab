%plots pdf overview with rf-summary-classification
%greschner
%2015-02-09 ath

clear all
% read in already done datasets

%rat
myfiles=dir('/Volumes/Lab/Projects/survey/rat/*.pdf');
cnt=1;
clear piecesRat
for i=1:length(myfiles)
    tmp = regexp(myfiles(i).name,'(\d+)-(\w+)-(\d+)','match');
    if cnt==1 || ~strcmp(piecesRat(cnt-1,:),tmp)
        piecesRat(cnt,:)=tmp;
        cnt=cnt+1;
    end
end

%mouse
myfiles=dir('/Volumes/Lab/Projects/survey/mouse/');
cnt=1;
clear piecesMouse
for i=1:length(myfiles)
    tmp = regexp(myfiles(i).name,'(\d+)-(\w+)-(\d+)','match');
    if ~isempty(tmp)
        if cnt==1 || ~strcmp(piecesMouse(cnt-1,:),tmp)
            piecesMouse(cnt,:)=tmp;
            cnt=cnt+1;
        end
    end
end


%monkey
file_path='/Volumes/Lab/Projects/survey/monkey/';
myfiles=dir('/Volumes/Lab/Projects/survey/monkey/*.pdf');
cnt=1;
clear piecesMonkey
for i=1:length(myfiles)
    tmp = regexp(myfiles(i).name,'(\d+)-(\w+)-(\d+)-(\d+)','match');
    if isempty(tmp)
        tmp = regexp(myfiles(i).name,'(\d+)-(\w+)-(\d+)','match');
    end
    if ~isempty(tmp)
        piecesMonkey(cnt,:)=tmp;
        cnt=cnt+1;
    end
end


% notebooks directory
notes=dir('/Volumes/Lab/Notebooks/');
expression='(\d+)-(\w+)-(\d+)';
clear mynote_date_list final_notes
cnt=1;
for i=1:length(notes)
    name=notes(i).name;
    if name(1)~='.' && name(1)~='#'
        mydate = regexpi(name,expression,'match');
        if ~isempty(mydate)
            fl=1;
            for j=1:length(piecesMouse)
                if strcmp(mydate,piecesMouse(j))
                    fl=0;
                    break;
                end
            end
            if fl
                for j=1:length(piecesRat)
                    if strcmp(mydate,piecesRat(j))
                        fl=0;
                        break;
                    end
                end
            end
            
            if fl
                myNotebook=['/Volumes/Lab/Notebooks/',name];
                if myNotebook(end)=='d'
                    myNotebook=[myNotebook,'/TXT.rtf'];
                end
                % attempt to read the notebook
                a=textread(myNotebook,'%s');
                
                for k=1:length(a)
                    if strcmpi(a{k},'rat') || strcmpi(a{k},'rat\') || strcmpi(a{k},'rat,') ...
                            || ~isempty(regexpi(a{k},'mouse')) || ~isempty(regexpi(a{k},'pig')) ...
                            || ~isempty(regexp(a{k},'squirrel'))
                        fl=0;
                        break;
                    end
                end
                if fl
                    mynote_date_list(cnt,:)=mydate;
                    final_notes{cnt}=name;
                    cnt=cnt+1;
                end
            end
        end
    end
end


% LOOK for modified data run folders, like in '2014-09-10-2/data000_nps'

analyzed=dir('/Volumes/Analysis/');
cnt=1;
clear myFiles
for i=1:length(analyzed)
    name=analyzed(i).name;
    if name(1)~='.'
        mydate = regexpi(name,'(\d+)-(\w+)-(\d+)','match');
        if ~isempty(mydate) % if it is a valid data folder
            for j=1:length(mynote_date_list)
                if strcmp(mydate{1},mynote_date_list{j}) % if notebook exists
                    
                    % check if 'data' exists
                    subfolders=dir(['/Volumes/Analysis/',name,'/data*']);
                    if isempty(subfolders)
                        subfolders=dir(['/Volumes/Analysis/',name,'/streamed/*']);
                    end
                    if ~isempty(subfolders)
                        myFiles{cnt,1}=name;
                        myFiles{cnt,2}=j; % number in final_notes
                        cnt=cnt+1;
                    end
                end
            end
        end
    end
end

expression='(\d+)-(\w+)-(\d+)-(\d+)';
secondExpr='(\d+)-(\w+)-(\d+)';
cnt=0;cnt1=1;
mycodes=cell(size(myFiles,1),4);
for i=1:size(myFiles,1)
    name=myFiles{i,1}
    mycodes{i,1}=name;
    
    mydate = regexpi(name,expression,'match');
    if isempty(mydate)
        mydate = regexpi(name,secondExpr,'match');
    end
    
    mydate=mydate{1};
    
    for j=1:length(piecesMonkey)
        if ~isempty(regexp(mydate,piecesMonkey{j}))
            mycodes{i,2}='done';
            break;
        end
    end
    
    if isempty(mycodes{i,2}) % not done yet
        
        myStim=[];
        
        % read notebook
        myNotebook=['/Volumes/Lab/Notebooks/',final_notes{myFiles{i,2}}];
        if myNotebook(end)=='d'
            myNotebook=[myNotebook,'/TXT.rtf'];
        end
        a=textread(myNotebook,'%s');
        
        
        % check subfolders for 'data0xx'
        subfolders=dir(['/Volumes/Analysis/',name,'/data*']);
        allruns={}; cnt=1;
        for j=1:length(subfolders)
            nameSF=subfolders(j).name;
            if length(nameSF)==7
                %check if there is dataxxx.params file
                tmp=dir(['/Volumes/Analysis/',name,'/',nameSF,'/',nameSF,'.params']);
                if ~isempty(tmp)
                    allruns{cnt}=nameSF;
                    cnt=cnt+1;
                    runtype='full';
                end
            end
        end
        
        %check streamed data if no analyzed runs found
        if isempty(allruns)
            subfolders=dir(['/Volumes/Analysis/',name,'/streamed/data*']);
            for j=1:length(subfolders)
                nameSF=subfolders(j).name;
                if length(nameSF)==7
                    %check if there is dataxxx.params file
                    tmp=dir(['/Volumes/Analysis/',name,'/streamed/',nameSF,'/',nameSF,'.params']);
                    if ~isempty(tmp)
                        allruns{cnt}=nameSF;
                        cnt=cnt+1;
                        runtype='streamed';
                    end
                end
            end
        end
        
        % check if found subfolders
        if isempty(allruns)
            mycodes{i,2}='No runs found';            
        else
            mycodes{i,3}=runtype;
            % find datarun in notebook
            beg=0;stop=0;
            for k=1:length(a)
                if ~beg && ~isempty(regexp(a{k},mydate, 'once'))
                    beg=k;
                elseif ~beg && ~isempty(regexp(a{k},[mydate,'.'], 'once'))
                    beg=k;
                elseif beg && isempty(regexp(a{k},mydate, 'once')) && ~isempty(regexp(a{k},[mydate(1:end-3),'.'], 'once'))
                    stop = k;
                    break;
                end
            end
            if stop<beg
                stop=length(a);
            end
            if beg==0 % if notebook was written for one piece, e.g. '2014-11-20-0.txt'
                beg=1;stop=length(a);
            end
            
            % parse notebook
            baseExpr='data.|Data.';
            WNexpr='BW.|RGB.|bw.|rgb.';
            myRun=[];myStim={};cntStim=1;
            for kk=beg+1:stop
                tmp=a{kk};
                
                tt=regexp(tmp,baseExpr); % check if it's data
                if ~isempty(tt) % if it is, check if it is one of the existing data runs
                    if ~isempty(regexp(tmp,'datarun\d.'))
                        tmp=[tmp(1:4),tmp(8:10)];
                    end
                    for ps=1:length(allruns)
                        tt1=regexp(tmp,[allruns{ps},'|','D',allruns{ps}(2:end)]);
                        if ~isempty(tt1)
                            myRun=allruns{ps};
                            break;
                        end
                        if ps==length(allruns)
                            myRun=[];
                        end
                    end
                end
                
                tt1=regexp(tmp,WNexpr); % check if it's RGB or BW run
                if isempty(tt1) % check if it is 'RGB bin xx' format
                    tt1=regexp(tmp,'BW|RGB|bw|rgb');
                    if ~isempty(tt1)
                        tt2=regexp(a{kk+1},'bin');
                        tt3=regexp(a{kk+2},'\d.','once');
                        tt4=regexp(a{kk+1},'\d.','once');
                        if ~isempty(tt2) && ~isempty(tt3)
                            stimCode=[a{kk},'-', a{kk+2}];
                        elseif ~isempty(tt4)
                            stimCode=[a{kk},'-', a{kk+1}];
                        end
                    end
                else
                    stimCode=a{kk};
                end
                if ~isempty(tt1) && ~isempty(myRun) && ~isempty(regexp(stimCode,'\d', 'once'))
                    myStim{cntStim,1}=myRun;
                    myStim{cntStim,2}=stimCode;
                    myRun=[];
                    cntStim=cntStim+1;
                end
            end
        end
        
        if ~isempty(myStim)
            mycodes{i,2}=myStim;
        else
            mycodes{i,2}='No WN runs found';
        end
    end
end




%% Create list

fid=fopen('/Volumes/Lab/Projects/survey/ath/list.txt','w');
for i=1:size(myFiles,1)
    
    if size(mycodes{i,2},2)==2 || size(mycodes{i,2},1)>1
        cols=[];sizes=[];
        for j=1:size(mycodes{i, 2},1)
            tmp=mycodes{i,2}{j,2};
            if ~isempty(regexpi(tmp,'bw','match'))
                cols(j)=1;
            elseif ~isempty(regexpi(tmp,'rgb','match'))
                cols(j)=3;
            end
            sizes(j)=str2num(regexpi(tmp,'\d+','match','once'));
            
        end
        
        if nnz(cols==3) && nnz(sizes(cols==3)>=6) % use RGB
            maxsize=max(sizes(cols==3));
            mydata=find(cols==3&sizes==maxsize,1);
            
        elseif nnz(cols==1) % use BW
            maxsize=max(sizes(cols==1));
            mydata=find(cols==1&sizes==maxsize,1);
        else % only small RGB runs
            maxsize=max(sizes);
            mydata=find(sizes==maxsize,1);
        end
        if strcmp(mycodes{i, 3},'streamed')
            s='/streamed/';
        else
            s='/';
        end
        path2data=['Volumes/Analysis/',mycodes{i,1},s,mycodes{i, 2}{mydata, 1}];
        fprintf(fid,'%s\r\n',path2data);
    end
end
fclose(fid);

%% plot something

file_path='/Volumes/Lab/Projects/survey/ath/';
file_name='list.txt';
ext_name='';

% open and load file
fid=fopen([file_path file_name]);
file=textscan(fid, '%s');
fclose(fid);
pathfile=file{1};

for i=1:length(pathfile)
    clear datarun
    t=dir(['/' pathfile{i} '/d*.params']);
    piece=regexp(pathfile{i},'(\d+)-(\w+)-(\d+)-(\d+)','match');
    if isempty(piece)
        piece=regexp(pathfile{i},'(\d+)-(\w+)-(\d+)','match');
    end
    
    for codeCNT=1:size(mycodes,1)
        if strcmp(mycodes{codeCNT,1},piece)
            codeN=codeCNT;
            break;
        end
    end
    
    if isempty(t)
        warning(['/' pathfile{i} '/*.params not found'])
    else
        datarun.names.rrs_params_path=['/' pathfile{i} '/' t.name];
        
        opt=struct('verbose',1,'load_params',1,'load_neurons',1);
        try datarun=load_data(datarun,opt);
            fl=1;
        catch err
            mycodes{codeN,4}=err;
            fl=0;
        end
        if fl
            info (datarun)
            
            paramsFile = edu.ucsc.neurobiology.vision.io.ParametersFile(datarun.names.rrs_params_path);
            
            
            class_list=[]; length_list=[];
            for ii=1:length(datarun.cell_types)
                if ~isempty(datarun.cell_types{ii}.name)
                    if isempty(regexpi(datarun.cell_types{ii}.name,'dupl'))...
                            &isempty(regexpi(datarun.cell_types{ii}.name,'junk'))...
                            &isempty(regexpi(datarun.cell_types{ii}.name,'crap'))...
                            &isempty(regexpi(datarun.cell_types{ii}.name,'weak'))...
                            &isempty(regexpi(datarun.cell_types{ii}.name,'contaminated'))...
                            &isempty(regexpi(datarun.cell_types{ii}.name,'edge'))...
                            &isempty(regexpi(datarun.cell_types{ii}.name,'bad'))...
                            &isempty(regexpi(datarun.cell_types{ii}.name,'unclassified'))...
                            &isempty(regexpi(datarun.cell_types{ii}.name,'nothing'))
                        class_list=[class_list ii];
                        length_list=[length_list length(datarun.cell_types{ii}.cell_ids)];
                    end
                end
            end
            
            if isempty(class_list)
                mycodes{codeN,4}='no cell classes'; % no cell classes
            elseif length(class_list)<3 && max(length_list)<10
                mycodes{codeN,4}='too few cells and classes';   % not enought cells in biggest class
            else
                mycodes{codeN,4}=class_list;
                
                ind=get_cell_indices(datarun,num2cell(class_list));
                t=zeros(length(ind),2);
                
                if isfield(datarun.vision.sta_fits{ind(1)},'mean')
                    for ii=1:length(ind)
                        t(ii,:)=datarun.vision.sta_fits{ind(ii)}.mean;
                    end
                    tt=sort(t(:,1));
                    x_range=[tt(2)-(tt(end-1)-tt(2))*.1 tt(end-1)+(tt(end-1)+tt(2))*.1];
                    %x_range=[min(t(:,1)) max(t(:,1))];
                    tt=sort(t(:,2));
                    y_range=[tt(2)-(tt(end-1)-tt(2))*.1 tt(end-1)+(tt(end-1)+tt(2))*.1];
                    %y_range=[min(t(:,2)) max(t(:,2))];
                    
                    fig=figure('Visible','off','PaperPosition',[0 0 12 8],'PaperSize',[12 8]);
                    set(gcf,'color','white');
                    
                    for ii=1:length(class_list)
                        
                        index=get_cell_indices(datarun,{class_list(ii)});
                        cell_id=datarun.cell_ids(index);
                        
                        if ~isempty(cell_id) && length(paramsFile.getArrayCell(cell_id(1),'GreenTimeCourse'))>0
                            ttime=zeros(length(cell_id),3,length(paramsFile.getArrayCell(cell_id(1),'GreenTimeCourse')));
                            tauto=zeros(length(cell_id),length(paramsFile.getArrayCell(cell_id(1),'Auto')));
                            for iii=1:length(cell_id)
                                t=paramsFile.getArrayCell(cell_id(iii), 'RedTimeCourse');
                                if ~isempty(t)
                                    ttime(iii,1,:)=t;
                                    ttime(iii,2,:)=paramsFile.getArrayCell(cell_id(iii), 'GreenTimeCourse');
                                    ttime(iii,3,:)=paramsFile.getArrayCell(cell_id(iii), 'BlueTimeCourse');
                                    tauto(iii,:)=paramsFile.getArrayCell(cell_id(iii), 'Auto');
                                end
                            end
                            if length(cell_id)>1
                                for iiii=1:3
                                    t=sum(abs(mean(ttime(:,iiii,:))));
                                    for iii=1:length(cell_id)
                                        ttime(iii,iiii,:)=ttime(iii,iiii,:)/sum(abs(ttime(iii,iiii,:)))*t;
                                    end
                                end
                                t=sum(abs(mean(tauto)));
                                for iii=1:length(cell_id)
                                    tauto(iii,:)=tauto(iii,:)/sum(abs(tauto(iii,:)))*t;
                                end
                            end
                            
                            
                            subplot(ceil(length(class_list)/2),6,(ii-1)*3+1)
                            plot_rf_fit(datarun, {class_list(ii)},'fill',false,'edge_color',[0 0 0]);
                            set(gca,'DataAspectRatio',[1 1 1]);
                            %set(gca,'PlotBoxAspectRatio',[1 1 1]);
                            %set(gca,'XTick',[],'YTick',[],'Visible','off');
                            set(gca,'XTick',[],'YTick',[]);
                            set(gca,'xLim',x_range,'yLim',y_range);
                            title(datarun.cell_types{class_list(ii)}.name);
                            %ylabel(datarun.cell_types{class_list(ii)}.name);
                            box on
                            subplot(ceil(length(class_list)/2),6,(ii-1)*3+2)
                            hold on
                            plot(squeeze(ttime(:,1,:))','r')
                            plot(squeeze(ttime(:,3,:))','b')
                            plot(squeeze(ttime(:,2,:))','g')
                            set(gca,'XTick',[],'YTick',[]);
                            set(gca,'xLim',[0  size(ttime,3)]);
                            box on;
                            subplot(ceil(length(class_list)/2),6,(ii-1)*3+3)
                            plot(tauto','k')
                            set(gca,'XTick',[],'YTick',[]);
                            set(gca,'xLim',[0  size(tauto,2)]);
                            box on;
                        end
                    end
                    
                    t=datarun.names.short_name;
                    if t(1)=='_', t=t(2:end); end
                    print(fig,'-dpdf',sprintf('%s%s%s.pdf',file_path,t,ext_name));
                    close(fig)
                else
                    mycodes{codeN,5}='no sta fits';
                end
            end
        end
        
    end
end








%% Martin's code

if 1
    file_name='dataset_overview.txt';
    ext_name='';
else
    file_name='dataset_overview_preliminary.txt';
    ext_name='_preliminary';
end

% open and load file
fid=fopen([file_path file_name]);
file=textscan(fid, '%s');
fclose(fid);
pathfile=file{1};

for i=22:length(pathfile)
    display(sprintf('%d - %d',i,length(pathfile)))
    
    clear datarun
    t=dir(['/' pathfile{i} '/d*.params']);
    if isempty(t)
        warning(['/' pathfile{i} '/*.params not found'])
    else
        datarun.names.rrs_params_path=['/' pathfile{i} '/' t.name];
        t=dir(['/' pathfile{i} '/rf-summary-classification.txt']);
        if ~isempty(t)
            datarun.names.rrs_classification_path=['/' pathfile{i} '/' t.name];
        end
        
        opt=struct('verbose',1,'load_params',1,'load_neurons',1);
        datarun=load_data(datarun,opt);
        datarun=load_vision_cell_types(datarun);
        info (datarun)
        
        paramsFile = edu.ucsc.neurobiology.vision.io.ParametersFile(datarun.names.rrs_params_path);
        
        fig=figure('Visible','off','PaperPosition',[0 0 12 8],'PaperSize',[12 8])
        set(gcf,'color','white');
        
        class_list=[1:4];
        for ii=5:length(datarun.cell_types)
            if ~isempty(datarun.cell_types{ii}.name)
                if isempty(strfind(datarun.cell_types{ii}.name,'duplicate'))&isempty(strfind(datarun.cell_types{ii}.name,'crap'))&isempty(strfind(datarun.cell_types{ii}.name,'weak'))&isempty(strfind(datarun.cell_types{ii}.name,'contaminated'))&isempty(strfind(datarun.cell_types{ii}.name,'edge'))
                    class_list=[class_list ii];
                end
            end
        end
        
        ind=get_cell_indices(datarun,num2cell(class_list));
        t=zeros(length(ind),2);
        for ii=1:length(ind)
            t(ii,:)=datarun.vision.sta_fits{ind(ii)}.mean;
        end
        tt=sort(t(:,1));
        x_range=[tt(2)-(tt(end-1)-tt(2))*.1 tt(end-1)+(tt(end-1)+tt(2))*.1];
        %x_range=[min(t(:,1)) max(t(:,1))];
        tt=sort(t(:,2));
        y_range=[tt(2)-(tt(end-1)-tt(2))*.1 tt(end-1)+(tt(end-1)+tt(2))*.1];
        %y_range=[min(t(:,2)) max(t(:,2))];
        
        for ii=1:length(class_list)
            
            index=get_cell_indices(datarun,{class_list(ii)});
            cell_id=datarun.cell_ids(index);
            
            if ~isempty(cell_id)
                ttime=zeros(length(cell_id),3,length(paramsFile.getArrayCell(cell_id(1),'GreenTimeCourse')));
                tauto=zeros(length(cell_id),length(paramsFile.getArrayCell(cell_id(1),'Auto')));
                for iii=1:length(cell_id)
                    t=paramsFile.getArrayCell(cell_id(iii), 'RedTimeCourse');
                    if ~isempty(t)
                        ttime(iii,1,:)=t;
                        ttime(iii,2,:)=paramsFile.getArrayCell(cell_id(iii), 'GreenTimeCourse');
                        ttime(iii,3,:)=paramsFile.getArrayCell(cell_id(iii), 'BlueTimeCourse');
                        tauto(iii,:)=paramsFile.getArrayCell(cell_id(iii), 'Auto');
                    end
                end
                if length(cell_id)>1
                    for iiii=1:3
                        t=sum(abs(mean(ttime(:,iiii,:))));
                        for iii=1:length(cell_id)
                            ttime(iii,iiii,:)=ttime(iii,iiii,:)/sum(abs(ttime(iii,iiii,:)))*t;
                        end
                    end
                    t=sum(abs(mean(tauto)));
                    for iii=1:length(cell_id)
                        tauto(iii,:)=tauto(iii,:)/sum(abs(tauto(iii,:)))*t;
                    end
                end
                
                
                subplot(ceil(length(class_list)/2),6,(ii-1)*3+1)
                plot_rf_fit(datarun, {class_list(ii)},'fill',false,'edge_color',[0 0 0]);
                set(gca,'DataAspectRatio',[1 1 1]);
                %set(gca,'PlotBoxAspectRatio',[1 1 1]);
                %set(gca,'XTick',[],'YTick',[],'Visible','off');
                set(gca,'XTick',[],'YTick',[]);
                set(gca,'xLim',x_range,'yLim',y_range);
                title(datarun.cell_types{class_list(ii)}.name);
                %ylabel(datarun.cell_types{class_list(ii)}.name);
                box on
                subplot(ceil(length(class_list)/2),6,(ii-1)*3+2)
                hold on
                plot(squeeze(ttime(:,1,:))','r')
                plot(squeeze(ttime(:,3,:))','b')
                plot(squeeze(ttime(:,2,:))','g')
                set(gca,'XTick',[],'YTick',[]);
                set(gca,'xLim',[0  size(ttime,3)]);
                box on;
                subplot(ceil(length(class_list)/2),6,(ii-1)*3+3)
                plot(tauto','k')
                set(gca,'XTick',[],'YTick',[]);
                set(gca,'xLim',[0  size(tauto,2)]);
                box on;
            end
        end
        
        t=datarun.names.short_name;
        t=strrep(t,'_marte','');
        t=strrep(t,'_jacob','');
        t=strrep(t,'__braid_snle_analysis-archive_Experiments_Array_Analysis_','');
        if t(1)=='_', t=t(2:end); end
        print(fig,'-dpdf',sprintf('%s%s%s.pdf',file_path,t,ext_name));
        close(fig)
    end
end





if 0
    fig = figure('Visible','off',...
        'PaperPosition',[0 0 6 4],...
        'PaperSize',[6 4])
    plot(rand(5))
    set(gca,'Position',[0 0 1 1])
    print(fig,'-dpdf',sprintf('%s%s.pdf',file_path,datarun.names.short_name))
    
    
    % first, create a multipage ps
    fnam='foo.ps';
    delete(fnam);
    for i=1:10
        plot(rand(10,1));
        print('-dpsc2','-append',fnam);
    end
end