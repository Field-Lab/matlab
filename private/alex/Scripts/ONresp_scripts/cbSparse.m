clear
date='20130703a'
% date='20130703b'

codeWord='cbSparse';

mainpath=['/Users/alexth/Desktop/old_stuff/',date,'/'];
hekapath=['/Users/alexth/Desktop/old_stuff/',date,'/HEKA/'];
heka=dir([hekapath,'*.phys']);

path2save=['/Users/alexth/Desktop/old_stuff/ONresp/',date,'/'];
if ~exist(path2save,'dir')
    mkdir(path2save);
end

% select heka files
clear ndFullList
file_list=[];
cnt=1;
for i=1:length(heka)
    if ~isempty(regexp(heka(i).name,codeWord, 'once'))
        file_list=[file_list i];
        
        a=regexp(heka(i).name,'ND');
        a=heka(i).name(a+2);
        if length(a)>1
            ndFullList(cnt)=str2num(a(1))+str2num(a(2));
        else
            ndFullList(cnt)=str2num(a(1));
        end
        cnt=cnt+1;
    end
end


% get protocol times
protocol=read_header_field_heka(hekapath, heka(file_list(1)).name, 'Stimulus Protocol');
protocol=round(protocol(2:end-1,[1,2]));

units=dir([mainpath,'units/*.mat']);
cbSparse=zeros(4500,length(file_list)*8,length(units));

names=cell(length(units),1);
for cnt=1:length(units)
    
    load([mainpath,'units/',units(cnt).name]);
    trialNumber=1;
    
    for i=1:length(file_list)
        spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
        conv=convolved(spikes,40,41500);
        conv=conv(121:end-120);
        for j=1:2:size(protocol,1)-1
            cbSparse(:,trialNumber,cnt)=conv(protocol(j,1)-499:protocol(j,1)+4000);
            trialNumber=trialNumber+1;
        end
        
    end  
    
    if isempty(regexp(units(cnt).name,date, 'once'))
        names{cnt}=[units(cnt).name(1),date,units(cnt).name(10:end-4)];
    else
        names{cnt}=units(cnt).name(1:end-4);
    end
end

save([path2save,'cbSparse'],'ndFullList','cbSparse','names')


%% Plot

clear

date='20130703b'
path2load=['/Users/alexth/Desktop/old_stuff/ONresp/',date,'/'];
load([path2load,'cbSparse'])
load([path2load,'quick'])
load('/Users/alexth/Desktop/old_stuff/_quickControl_/CB_sparse.mat')


coord_pos=cell(10,10);
coord_neg=cell(10,10);
for i=1:10
    for j=1:10
        coord_pos{i,j}=find(matr(i,j,:)==1);
        coord_neg{i,j}=find(matr(i,j,:)==-1);
    end
end

path2save=['/Users/alexth/Desktop/old_stuff/ONresp/',date,'/cbSparse_plot/'];
if ~exist(path2save,'dir')
    mkdir(path2save);
end

figure
set(gcf,'position',[4          90        1259         859])
for k=1:size(cbSparse,3)

    m=1;
    for i=1:10
        for j=1:10
            subplot(10,10,m)
            a=mean(cbSparse(1:3500,coord_neg{i,j},k),2);
            b=mean(cbSparse(1:3500,coord_pos{i,j},k),2);
            hold off
            plot(a)
            hold on
            plot(b,'r')
            
            a=mean(cbSparse(1:3500,coord_neg{i,j}+152,k),2);
            b=mean(cbSparse(1:3500,coord_pos{i,j}+152,k),2);
            plot(a,'g')
            hold on
            plot(b,'m')
            
            line([500 500],[0 35],'color','k')
            line([2500 2500],[0 35],'color','k')
            axis([1 3500 0 35])
            set(gca,'xtick',0,'xticklabel','')
            m=m+1;
        end
    end
    
    subplot('position',[0.92 0.87 0.07 0.07])
    hold off
    plot(mean(black_flash(:,17:26,k),2),'b')
    hold on
    plot(mean(white_flash(:,17:26,k),2),'r')
    
    plot(mean(black_flash(:,35:44,k),2),'g')
    hold on
    plot(mean(white_flash(:,35:44,k),2),'m')
    
    line([500 500],[0 50],'color','k')
    line([2500 2500],[0 50],'color','k')
    axis([1 3500 0 50])

    
    subplot('Position',[0.5 0.95 0.00001 0.00001])
    set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
    title([names{k},'  Black, White'],'FontSize',12,'FontWeight','bold','Interpreter','None')
    saveas(gcf,[path2save,names{k},'.png'])
end




%% Plot for Thomas

clear

date='20130703a'
path2load=['/Users/alexth/Desktop/old_stuff/ONresp/',date,'/'];
load([path2load,'cbSparse'])
load([path2load,'quick'])
load('/Users/alexth/Desktop/old_stuff/_quickControl_/CB_sparse.mat')


coord_pos=cell(10,10);
coord_neg=cell(10,10);
for i=1:10
    for j=1:10
        coord_pos{i,j}=find(matr(i,j,:)==1);
        coord_neg{i,j}=find(matr(i,j,:)==-1);
    end
end

% path2save=['/Users/alexth/Desktop/old_stuff/ONresp/',date,'/cbSparse_plot_for_Thomas/'];
path2save=['/Users/alexth/Desktop/old_stuff/ONresp/cbSparse_plot_for_Thomas/'];
if ~exist(path2save,'dir')
    mkdir(path2save);
end

%%%%%% Plot Both %%%%%%

for k=1:size(cbSparse,3)

    %white flash
    list_pos_w=zeros(5,4);
    m=1;
    cnt=1;
    for i=1:10
        for j=1:10
%             subplot(10,10,m)
            data=cbSparse(1:3500,coord_pos{i,j},k);
            a=corr(data);
            % checks if there are only few NAN entries and if the rest is
            % well correlated
            if sum(sum(triu(a)>0.5&triu(a)<1))>(size(a,2)-1)          
                list_pos_w(cnt,1:2)=[i,j];
                cnt=cnt+1;
            end
            
            
            data=cbSparse(1:3500,coord_pos{i,j}+152,k);
            a=corr(data);
            % checks if there are only few NAN entries and if the rest is
            % well correlated
            if sum(sum(triu(a)>0.5&triu(a)<1))>(size(a,2)-1)          
                 list_pos_w(cnt,3:4)=[i,j];
            end
            m=m+1;
        end
    end
   
     %black flash
    list_pos_b=zeros(5,4);
    m=1;
    cnt=1;
    for i=1:10
        for j=1:10
%             subplot(10,10,m)
            data=cbSparse(1:3500,coord_neg{i,j},k);
            a=corr(data);
            % checks if there are only few NAN entries and if the rest is
            % well correlated
            if sum(sum(triu(a)>0.5&triu(a)<1))>(size(a,2)-1)          
                list_pos_b(cnt,1:2)=[i,j];
                cnt=cnt+1;
            end
            
            
            data=cbSparse(1:3500,coord_neg{i,j}+152,k);
            a=corr(data);
            % checks if there are only few NAN entries and if the rest is
            % well correlated
            if sum(sum(triu(a)>0.5&triu(a)<1))>(size(a,2)-1)          
                 list_pos_b(cnt,3:4)=[i,j];
            end
            m=m+1;
        end
    end
 
    if sum(list_pos_b(:))||sum(list_pos_w(:))
        figure
        set(gcf,'position',[4   195   634   735])
        cnt=1;
        for m=1:5
            if list_pos_w(m,1)
                subplot(5,2,cnt)
                i=list_pos_w(m,1);j=list_pos_w(m,2);
                plot(mean(cbSparse(1:3500,coord_pos{i,j},k),2),'m','linewidth',2)
                hold on
                plot(3701:7200,mean(cbSparse(1:3500,coord_neg{i,j},k),2),'c','linewidth',2)
                
                plot(mean(white_flash(1:3500,17:26,k),2),'r')
                plot(3701:7200,mean(black_flash(1:3500,17:26,k),2),'b')
                
                axis tight
                p=get(gca,'YLim');
                line([500 500],[0 p(2)],'color','k')
                line([2500 2500],[0 p(2)],'color','k')
                line([500 500]+3701,[0 p(2)],'color','k')
                line([2500 2500]+3701,[0 p(2)],'color','k')
                title(['ND6, position i=',int2str(i),', j=',int2str(j)])
            end
            
            if list_pos_w(m,3)
                subplot(5,2,cnt+1)
                i=list_pos_w(m,3);j=list_pos_w(m,4);
                plot(mean(cbSparse(1:3500,coord_pos{i,j}+152,k),2),'m','linewidth',2)
                hold on
                plot(3701:7200,mean(cbSparse(1:3500,coord_neg{i,j}+152,k),2),'c','linewidth',2)
                
                plot(mean(white_flash(1:3500,35:44,k),2),'r')
                plot(3701:7200,mean(black_flash(1:3500,35:44,k),2),'b')
                
                axis tight
                p=get(gca,'YLim');
                line([500 500],[0 p(2)],'color','k')
                line([2500 2500],[0 p(2)],'color','k')
                line([500 500]+3701,[0 p(2)],'color','k')
                line([2500 2500]+3701,[0 p(2)],'color','k')
                title(['ND4, position i=',int2str(i),', j=',int2str(j)])
            end
            cnt=cnt+2;
            
        end
        
        subplot('Position',[0.5 0.95 0.00001 0.00001])
        set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
        title([names{k},'  White, Black'],'FontSize',12,'FontWeight','bold','Interpreter','None')
        saveas(gcf,[path2save,names{k},'.png'])
        close all
    end
end


%%%%% ND6 %%%%%

for k=1:size(cbSparse,3)

    %white flash
    rightspot=zeros(10);
    list_pos_w=[];
    m=1;
    for i=1:10
        for j=1:10
%             subplot(10,10,m)
            data=cbSparse(1:3500,coord_pos{i,j},k);
            a=corr(data);
            % checks if there are only few NAN entries and if the rest is
            % well correlated
            if sum(sum(triu(a)>0.5&triu(a)<1))>(size(a,2)-1)          
                rightspot(i,j)=1;
                list_pos_w=[list_pos_w; [i,j]];
%                 plot(data,'r')
            else
                rightspot(i,j)=0;
%                 plot(data,'b')
            end

%             axis tight           
%             line([500 500],[0 35],'color','k')
%             line([2500 2500],[0 35],'color','k')
%             set(gca,'xtick',0,'xticklabel','')
            m=m+1;
        end
    end
   
     %black flash
    rightspot=zeros(10);
    list_pos_b=[];
    m=1;
    for i=1:10
        for j=1:10
%             subplot(10,10,m)
            data=cbSparse(1:3500,coord_neg{i,j},k);
            a=corr(data);
            % checks if there are only few NAN entries and if the rest is
            % well correlated
            if sum(sum(triu(a)>0.5&triu(a)<1))>(size(a,2)-1)          
                rightspot(i,j)=1;
%                  plot(data,'m')
                list_pos_b=[list_pos_b; [i,j]];
            else
                rightspot(i,j)=0;
%                 plot(data,'c')
            end

%             axis tight           
%             line([500 500],[0 35],'color','k')
%             line([2500 2500],[0 35],'color','k')
%             set(gca,'xtick',0,'xticklabel','')
            m=m+1;
        end
    end
 
    if size(list_pos_w,1)
        figure
        set(gcf,'position',[4   195   634   735])
        for m=1:size(list_pos_w,1)
            
            i=list_pos_w(m,1);j=list_pos_w(m,2);
            subplot(size(list_pos_w,1),1,m)
%             ip=i-min(list_pos_w(m,1))+1;jp=j-min(list_pos_w(m,2))+1;
%             subplot(10,10,(i-1)*10+j)
            
            plot(mean(cbSparse(1:3500,coord_pos{i,j},k),2),'m','linewidth',2)
            hold on
            plot(3701:7200,mean(cbSparse(1:3500,coord_neg{i,j},k),2),'c','linewidth',2)
            
            plot(mean(white_flash(1:3500,17:26,k),2),'r')
            plot(3701:7200,mean(black_flash(1:3500,17:26,k),2),'b')
            
            axis tight
            p=get(gca,'YLim');
            line([500 500],[0 p(2)],'color','k')
            line([2500 2500],[0 p(2)],'color','k')
            line([500 500]+3701,[0 p(2)],'color','k')
            line([2500 2500]+3701,[0 p(2)],'color','k')
        end
        
        subplot('Position',[0.5 0.95 0.00001 0.00001])
        set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
        title([names{k},'  White, Black, ND6'],'FontSize',12,'FontWeight','bold','Interpreter','None')
        saveas(gcf,[path2save,names{k},'_ND6.png'])
        close all
    end
end


%%%%% ND4 %%%%%

for k=1:size(cbSparse,3)

    %white flash
    rightspot=zeros(10);
    list_pos_w=[];
    m=1;
    for i=1:10
        for j=1:10
%             subplot(10,10,m)
            data=cbSparse(1:3500,coord_pos{i,j}+152,k);
            a=corr(data);
            % checks if there are only few NAN entries and if the rest is
            % well correlated
            if sum(sum(triu(a)>0.5&triu(a)<1))>(size(a,2)-1)          
                rightspot(i,j)=1;
                list_pos_w=[list_pos_w; [i,j]];
%                 plot(data,'r')
            else
                rightspot(i,j)=0;
%                 plot(data,'b')
            end

%             axis tight           
%             line([500 500],[0 35],'color','k')
%             line([2500 2500],[0 35],'color','k')
%             set(gca,'xtick',0,'xticklabel','')
            m=m+1;
        end
    end
   
     %black flash
    rightspot=zeros(10);
    list_pos_b=[];
    m=1;
    for i=1:10
        for j=1:10
%             subplot(10,10,m)
            data=cbSparse(1:3500,coord_neg{i,j}+152,k);
            a=corr(data);
            % checks if there are only few NAN entries and if the rest is
            % well correlated
            if sum(sum(triu(a)>0.5&triu(a)<1))>(size(a,2)-1)          
                rightspot(i,j)=1;
%                 plot(data,'m')
                list_pos_b=[list_pos_b; [i,j]];
            else
                rightspot(i,j)=0;
%                 plot(data,'c')
            end

%             axis tight           
%             line([500 500],[0 35],'color','k')
%             line([2500 2500],[0 35],'color','k')
%             set(gca,'xtick',0,'xticklabel','')
            m=m+1;
        end
    end
 
    if size(list_pos_w,1)
        figure
        set(gcf,'position',[4   195   634   735])
        for m=1:size(list_pos_w,1)
            subplot(size(list_pos_w,1),1,m)
            i=list_pos_w(m,1);j=list_pos_w(m,2);
            plot(mean(cbSparse(1:3500,coord_pos{i,j}+152,k),2),'m','linewidth',2)
            hold on
            plot(3701:7200,mean(cbSparse(1:3500,coord_neg{i,j}+152,k),2),'c','linewidth',2)
            
            plot(mean(white_flash(1:3500,35:44,k),2),'r')
            plot(3701:7200,mean(black_flash(1:3500,35:44,k),2),'b')
            
            axis tight
            p=get(gca,'YLim');
            line([500 500],[0 p(2)],'color','k')
            line([2500 2500],[0 p(2)],'color','k')
            line([500 500]+3701,[0 p(2)],'color','k')
            line([2500 2500]+3701,[0 p(2)],'color','k')
        end
        
        subplot('Position',[0.5 0.95 0.00001 0.00001])
        set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
        title([names{k},'  White, Black, ND4'],'FontSize',12,'FontWeight','bold','Interpreter','None')
        saveas(gcf,[path2save,names{k},'_ND4.png'])
        close all
    end
end
