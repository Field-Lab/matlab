cd('/mnt/muench_data/user/alexandra/scripts')
clear

dates=cell(15,1);
dates{1}='20130220';
dates{2}='20130220_1';
dates{3}='20130224';
dates{4}='20130225';
dates{5}='20130226';
dates{6}='20130227';
dates{7}='20130301';
dates{8}='20130301_1';
dates{9}='20130301_2';
dates{10}='20130302';
dates{11}='20130302_1';
dates{12}='20120329';
dates{13}='20120627';
dates{14}='20120714';
dates{15}='20121023';

par=4;
hc=cell(1,8);
lc=cell(1,8);
for datesCNT=1:15
    date=dates{datesCNT}
    path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];    
    load([path2save,date,'_nonl'])
    load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date],'names','onOff')
    load([path2save,'gain'],'nonLinear','form')
    if length(unique(form))==2
        highC=nonLinear(par,1:2:end,:);
        lowC=nonLinear(par,2:2:end,:);
        
        a=reshape(highC,size(highC,2),size(highC,3));
        b=reshape(lowC,size(lowC,2),size(lowC,3));
        tr=size(highC,2)/8;
        cc=1;
        if datesCNT==14
            tr=size(highC,2)/7;
            br=nanmean(a(1:tr,:));
            hc{cc}=[hc{cc} zeros(size(br))];
            lc{cc}=[lc{cc} zeros(size(br))];
            cc=2;
        end
        for i=1:tr:size(highC,2)
            hc{cc}=[hc{cc} nanmean(a(i:i+tr-1,:))];
            lc{cc}=[lc{cc} nanmean(b(i:i+tr-1,:))];
            cc=cc+1;
        end
    else
        highC=nonLinear(par,:,:);
        a=reshape(highC,size(highC,2),size(highC,3));
        tr=size(highC,2)/8;
        cc=1;
        for i=1:tr:size(highC,2)
            b=nanmean(a(i:i+tr-1,:));
            hc{cc}=[hc{cc} b];
            lc{cc}=[lc{cc} zeros(size(b))];
            cc=cc+1;
        end
    end
end

hc_gain=zeros(520,8);
lc_gain=hc_gain;
for i=1:8
    hc_gain(1:520,i)=hc{i};
    lc_gain(1:520,i)=lc{i};
end

% save('/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_gain','hc_gain','lc_gain')
% load('/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_all','onOff','bothContr')


load('/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_all','exp_codes','bothContr')


%par1
for i=1:8
    b=hc_gain(exp_codes>5,i);
    b(b<=-10)=nan;
    b(b>40)=nan;
    par1(i,1)=nanmean(b);
    
    b=hc_gain(exp_codes<6,i);
    b(b<=-10)=nan;
    b(b>40)=nan;
    par1(i,2)=nanmean(b);
    
    b=lc_gain(exp_codes>5,i);
    b(b<=-10)=nan;
    b(b>40)=nan;
    par1(i,3)=nanmean(b);

end
%par2
for i=1:8
    b=hc_gain(exp_codes>5,i);
    b(b<=0)=nan;
    b(b>100)=nan;
    par2(i,1)=nanmean(b);
    
    b=hc_gain(exp_codes<6,i);
    b(b<=0)=nan;
    b(b>100)=nan;
    par2(i,2)=nanmean(b);
    
    b=lc_gain(exp_codes>5,i);
    b(b<=0)=nan;
    b(b>100)=nan;
    par2(i,3)=nanmean(b);

end
%par3
for i=1:8
    b=hc_gain(exp_codes>5,i);
    b(b<=0)=nan;
    b(b>0.2)=nan;
    par3(i,1)=nanmean(b);
    
    b=hc_gain(exp_codes<6,i);
    b(b<=0)=nan;
    b(b>0.2)=nan;
    par3(i,2)=nanmean(b);
    
    b=lc_gain(exp_codes>5,i);
    b(b<=0)=nan;
    b(b>0.2)=nan;
    par3(i,3)=nanmean(b);

end
%par4
for i=1:8
    b=hc_gain(exp_codes>5,i);
    b(b<=0)=nan;
    b(b>0.15)=nan;
    par4(i,1)=nanmean(b);
    
    b=hc_gain(exp_codes<6,i);
    b(b<=0)=nan;
    b(b>0.15)=nan;
    par4(i,2)=nanmean(b);
    
    b=lc_gain(exp_codes>5,i);
    b(b<=0)=nan;
    b(b>0.15)=nan;
    par4(i,3)=nanmean(b);

end
      
m=[0.07 0.025 0.02 0.05 0.055 0.05 0.07 0.07]+0.02;
gain_k=3;
baseline_k=0;
x=rangeHC;
figure
for i=1:8
    subplot(2,4,i)
    b=(1+exp((par3(i,1)-m(i)-x)/(par4(i,1)/gain_k)));
    y=par1(i,1)*baseline_k+ones(size(b)).*par2(i,3)./b;
    plot(x,y)
    hold on
end
% for i=1:8
%     subplot(2,4,i)
%     b=(1+exp((par3(i,2)-x)/par4(i,2)));
%     y=par1(i,2)+ones(size(b)).*par2(i,2)./b;
%     plot(x,y,'k')
% end
x=rangeHC;
for i=1:8
    subplot(2,4,i)
    b=(1+exp((par3(i,3)-x)/par4(i,3)));
    y=par1(i,3)*baseline_k+ones(size(b)).*par2(i,3)./b;
    plot(x,y,'r')
    axis([-0.4 0.4 0 50])
end





m=[0.07 0.025 0.02 0.05 0.055 0.05 0.07 0.07]+0.02;
m=zeros(8,1);
gain_k=1;
baseline_k=1;
x=rangeHC;
figure
for i=1:8
    subplot(2,4,i)
    b=(1+exp((par3(i,1)-m(i)-x)/(par4(i,1)/gain_k)));
    y=par1(i,1)*baseline_k+ones(size(b)).*par2(i,1)./b;
    plot(x,y)
    hold on
end

x=rangeLC;
for i=1:8
    subplot(2,4,i)
    b=(1+exp((par3(i,3)-x)/par4(i,3)));
    y=par1(i,3)*baseline_k+ones(size(b)).*par2(i,3)./b;
    plot(x,y,'r')
    axis([-0.4 0.4 0 50])
end


g = fittype('a+b/(1+exp((c-x)/d))');








for i=1:8
    subplot(2,4,i)
    a=hc_gain(:,i);
%     a(a==0)=nan;
    plot(a,'.','MarkerSize',5)
    hold on
    b=lc_gain(:,i);    
    plot(b,'.r','MarkerSize',5)
    
    axis([0 520 -0.02 0.12])
end

figure
for i=1:8
    subplot(2,4,i)
    a=lc_gain(onOff<0,i);
%     a(a==0)=nan;
    plot(a,'.')
    hold on
    b=lc_gain(onOff>0,i);    
    plot(b,'.r')
    
    axis([0 520 -0.02 0.12])
end


figure
for i=1:8
    subplot(2,4,i)
    a=hc_gain(onOff<0,i)./lc_gain(onOff<0,i);
%     a(a==0)=nan;
    plot(a,'.')
    hold on
    a=hc_gain(onOff>0,i)./lc_gain(onOff>0,i);
    plot(a,'.r')
    
    axis([0 249 0 15])
end





















cd('/mnt/muench_data/user/alexandra/scripts')
clear
date='20121023'
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date],'path2save')


rangeHC=-0.3:0.3/30:0.3;
rangeLC=-0.07:0.07/30:0.07;

load([path2save,date,'_nonlinear'],'nonLinear','names')
nonlin=zeros(min(size(LinearFilter,3),24*8),61,length(units));


for i=1:min(size(LinearFilter,3),24*8)
    i
    for t=1:dim
        subflicker=flicker(correctedProtocols(per*(t-1)+1,1,i)-startTimes(i)+1:correctedProtocols(per*(t-1)+dur,1,i)-startTimes(i),i);
        tmp=zeros(length(subflicker)-499,length(units));
        a=reshape(LinearFilter(500:-1:1,t,i,:),500,length(units))';
        for kk=1:size(a,1)
            a(kk,:)=a(kk,:)-mean(a(kk,:));
            a(kk,:)=a(kk,:)./sum(abs(a(kk,:)));
        end
        for rt=500:length(subflicker)
            tmp(rt-499,:)=a*subflicker(rt-499:rt);
        end
        
        a=convSpikes((correctedProtocols(per*(t-1)+1,1,i)-startTimes(i))+500:(correctedProtocols(per*(t-1)+dur,1,i)-startTimes(i)),i,:);
        a=reshape(a,size(a,1),size(a,3));
        
        if std(subflicker)>0.15
            x=rangeHC;
            stepSize=0.3/30;
        else
            x=rangeLC;
            stepSize=0.07/30;
        end
        
        kk=1;
        for rt=x
            tt=tmp>=rt&tmp<rt+stepSize;
            [timeP,cellID]=find(tt);
            for l=unique(cellID)'
                nonlin(i,kk,l)=mean(a(timeP(cellID==l),l));
            end
            kk=kk+1;
        end
    end
end
save([path2save,date,'_nonlin'],'nonlin','names')

nonl=zeros(64,61,length(units));
cc=1;
nonlin_t=nonlin;
nonlin_t(nonlin_t==0)=nan;
for i=1:6:min(size(LinearFilter,3),24*8)
    nonl(cc,:,:)=nanmean(nonlin_t(i:2:i+5,:,:));
    nonl(cc+1,:,:)=nanmean(nonlin_t(i+1:2:i+5,:,:));
    cc=cc+2;
end

g = fittype('a+b/(1+exp((c-x)/d))');
for i=1:64
    i
    if mod(i,2)==1
        x=rangeHC;
    else
        x=rangeLC;
    end
    for cnt=1:size(LinearFilter,4)
        a=nonl(i,:,cnt);
        b=find(isnan(a(1:30)));
        if ~isempty(b)
            a(b)=min(a(b(end)+1:b(end)+4));
        end
        b=find(isnan(a(30:end)))+29;
        if ~isempty(b)            
            a(b)=max(a(b(1)-8:b(1)-1));
        end
%         a=a(5:end-5);
%         x=x(5:end-5);

        try
            tmp_f=max(a)-min(a);
            tmp_f=find(a>tmp_f,1);
            res=fit(x',a',g,'StartPoint', [min(a),max(a),x(tmp_f),0.03]);
            nonLinear(1,i,cnt)=res.a;
            nonLinear(2,i,cnt)=res.b;
            nonLinear(3,i,cnt)=res.c;
            nonLinear(4,i,cnt)=res.d;
        catch
            disp('could NOT fit')
            date
            cnt
            t
            disp('END OF could NOT fit')
        end
    end
end  
         


figure
for i=10:2:16
    if mod(i,2)==1
        x=rangeHC;
    else
        x=rangeLC;
    end
    for cnt=1:39
        subplot(7,6,cnt)
        a=nonl(i,:,cnt);
        b=find(isnan(a(1:30)));
        if ~isempty(b)
            a(b)=min(a(b(end)+1:b(end)+4));
        end
        b=find(isnan(a(30:end)))+29;
        if ~isempty(b)
            a(b)=max(a(b(1)-8:b(1)-1));
        end
        plot(x,a)
        hold on
        y_t=nonLinear(1,i,cnt)+ones(size(x))*nonLinear(2,i,cnt)./(1+exp((nonLinear(3,i,cnt)-x)/nonLinear(4,i,cnt)));
        plot(x,y_t,'r')
        drawnow
    end
end











