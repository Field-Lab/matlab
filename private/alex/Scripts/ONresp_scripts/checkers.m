clear
cd('/Users/alexth/Desktop/scripts/ONresp_scripts')

% date='20130703a'
date='20130703b'

codeWord='check';

mainpath=['/Users/alexth/Desktop/old_stuff/',date,'/'];
hekapath=['/Users/alexth/Desktop/old_stuff/',date,'/HEKA/'];
heka=dir([hekapath,'*.phys']);


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

units=dir([mainpath,'units/*.mat']);
names=cell(length(units),1);

load('/Users/alexth/Desktop/old_stuff/_quickControl_/s.mat')
full_matr=zeros(18000,1600);
defaultStream = RandStream.getDefaultStream;
defaultStream.State = s(1).State;
for i=1:18000
    full_matr(i,:)=uint8(randi([0 1],1,1600));
end
full_matr=full_matr';
full_matr(full_matr<1)=-1;
clear b s defaultStream

path2save=['/Users/alexth/Desktop/old_stuff/ONresp/',date,'/checkers_plot/'];
if ~exist(path2save,'dir')
    mkdir(path2save);
end

figure
set(gcf,'position',[105         204        1145         727])

for cnt=1:length(units)
    
    load([mainpath,'units/',units(cnt).name]);
    
    if isempty(regexp(units(cnt).name,date, 'once'))
        names{cnt}=[units(cnt).name(1),date,units(cnt).name(10:end-4)];
    else
        names{cnt}=units(cnt).name(1:end-4);
    end
    
    
    for i=1
        protocol=read_header_field_heka(hekapath, heka(file_list(i)).name, 'Stimulus Protocol');
        protocol=round(protocol(1:end,1:4));
        load([mainpath,'checker/',date(1:8),'_#',heka(file_list(i)).name(13:16),'_checker_h_seq_.mat'])
        FlipTimeStamps=FlipTimeStamps-FlipTimeStamps(1);
        FlipTimeStamps=FlipTimeStamps*1000;
        
        spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
        spikes=spikes-protocol(3,1);
        spikes=round(spikes/(100/3));
        spikes(spikes<16)=[];
        tt=length(spikes);
        a=zeros(1,spikes(end));
        
        while sum(spikes)>0
            [~,k]=unique(spikes);
            m=spikes(k);
            m(m==0)=[];
            a(m)=a(m)+1;
            spikes(k)=0;            
        end        

        LF=zeros(1600,15);
        
        for m=1:max(unique(a))
            b=find(a==m);
            for j=1:15
                LF(:,j)=LF(:,j)+sum(full_matr(:,b-j+1),2);
            end
        end 
        LF=LF/tt;

        for j=1:15
            subplot(3,5,j)
            k=reshape(LF(:,j),40,40);
            k(1,1)=max(LF(:));
            k(2,1)=min(LF(:));
            imagesc(k);
            colormap gray            
        end
        subplot('Position',[0.5 0.95 0.00001 0.00001])
        set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
        title([names{cnt},'   ND6'],'FontSize',12,'FontWeight','bold','Interpreter','None')
        saveas(gcf,[path2save,names{cnt},'.png'])        
    end  
    

end

