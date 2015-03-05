cd('/mnt/muench_data/user/alexandra/scripts')
clear

dates=cell(6,1);
dates{1}='20130220';
dates{2}='20130220_1';
dates{3}='20130224';
dates{4}='20130225';
dates{5}='20130226';
dates{6}='20121023';
k=[];
for i=1:24:24*8
    k=[k i:2:i+17];
end

names_all=[];
onOff_all=[];
exp_codes=[];
nonl=zeros(61,72,300);
last=1;
for datesCNT=1:6
    date=dates{datesCNT}
    path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];
    load([path2save,date,'_nonlin'],'nonlin','names','rangeHC','rangeLC')
    load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date],'units','onOff')
    
    names_all=[names_all names];
    exp_codes=[exp_codes ones(1,size(nonlin,4))*datesCNT];
    onOff_all=[onOff_all onOff];
    
    
    if size(nonlin,3)==192
        nonlin=reshape(nonlin,61,192,size(nonlin,4));
        nonl(:,:,last:last+size(nonlin,3)-1)=nonlin(:,k,:);
    else
        strt=10;
        if strcmp(date,'20130220')
            strt=11;
        end
        nonlin=reshape(nonlin,61,size(nonlin,3),size(nonlin,4));
        nonl(:,:,last:last+size(nonlin,3)-1)=nonlin(:,strt:strt+71,:);
    end
    last=last+size(nonlin,3);
end
nonl=nonl(:,:,1:last-1);
names=names_all;
onOff=onOff_all;
save('/mnt/muench_data/data/alexandra/MEA_data/analysis/temporal_develop_nonl.mat','nonl','names','onOff','rangeHC','exp_codes')


clear
dates=cell(6,1);
dates{1}='20130220';
dates{2}='20130220_1';
dates{3}='20130224';
dates{4}='20130225';
dates{5}='20130226';
dates{6}='20121023';
g = fittype('a+b/(1+exp((c-x)/d))');
load('/mnt/muench_data/data/alexandra/MEA_data/analysis/temporal_develop_nonl.mat');
nonLinear=zeros(4,size(nonl,2),size(nonl,3));
x=rangeHC;
for i=1:size(nonl,2)
    i

    for cnt=1:size(nonl,3)
        a=nonl(:,i,cnt);
        if sum(isnan(a))<30
            b=find(isnan(a(1:30)));
            d=find(a(1:30)==0);
            if isempty(b)||b(end)<d(end)
                b=d;
            end
            if ~isempty(b)
                a(b)=min(a(b(end)+1:b(end)+2));
            end
            
            b=find(isnan(a(30:end)))+29;
            if ~isempty(b)
                c=sort(a(~isnan(a)));
                c=mean(c(end-5:end));
                a(b)=c;
            end
            d=find(a(30:end)==0)+29;
            if ~isempty(d)
                c=mean(a(d(1)-5:d(1)-1));
                 a(d)=c;
            end
            
            try
                tmp_f=max(a)-min(a);
                tmp_f=find(a>tmp_f,1);
                res=fit(x',a,g,'StartPoint', [min(a),max(a),x(tmp_f),0.03]);
                nonLinear(1,i,cnt)=res.a;
                nonLinear(2,i,cnt)=res.b;
                nonLinear(3,i,cnt)=res.c;
                nonLinear(4,i,cnt)=res.d;
            catch
                disp('could NOT fit')
                date
                cnt
                disp('END OF could NOT fit')
            end
        end
    end
end
save('/mnt/muench_data/data/alexandra/MEA_data/analysis/temporal_develop_gain','nonLinear','names','onOff','exp_codes')





for i=1:72
    tmp=nonLinear(4,i,:);
    tmp=ones(size(tmp))./tmp;
    tmp=tmp(tmp>0&tmp<20);
    k(i)=nanmean(tmp);
    m(i)=nanstd(tmp)/sqrt(length(tmp));
end

j=1;
for i=1:9:72
    errorbar(j:j+8,k(i:i+8),m(i:i+8))
    j=j+12;
    hold on
end
axis([0 98 0 20])
set(gca,'xtick',5:12:100,'xticklabel',{'8','7','6','5','4','3','2','1'})
xlabel('ND')
ylabel('gain, a.u.')
title('Fig.M4. Temporal evolution of gain at high contrast GWN')



