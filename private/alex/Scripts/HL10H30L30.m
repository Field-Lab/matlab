cd('/mnt/muench_data/user/alexandra/scripts')

clear
codeWord='HL10';
load(['/mnt/muench_data/data/alexandra/MEA_data/sine_analysis/',codeWord,'_accumParameters'])
a=who;
for i=1:size(a,1)
    if ~strcmp(a{i},'codeWord')
        eval([a{i},'_',codeWord,'=',a{i},';']);
    end
end
clear i
codeWord='H30s';
load(['/mnt/muench_data/data/alexandra/MEA_data/sine_analysis/',codeWord,'_accumParameters'])
a=who;
for i=1:size(a,1)
    if ~strcmp(a{i},'a')&&~strcmp(a{i},'codeWord')&&isempty(regexp(a{i},'HL10', 'once'))
        eval([a{i},'_',codeWord,'=',a{i},';']);
    end
end
clear i
codeWord='L30s';
load(['/mnt/muench_data/data/alexandra/MEA_data/sine_analysis/',codeWord,'_accumParameters'])
a=who;
for i=1:size(a,1)
    if ~strcmp(a{i},'a')&&~strcmp(a{i},'codeWord')&&isempty(regexp(a{i},'HL10', 'once'))&&isempty(regexp(a{i},'H30s', 'once'))
        eval([a{i},'_',codeWord,'=',a{i},';']);
    end
end

a=who;
for i=1:size(a,1)
    if ~strcmp(a{i},'a')&&isempty(regexp(a{i},'HL10', 'once'))&&isempty(regexp(a{i},'H30s', 'once'))&&isempty(regexp(a{i},'L30s', 'once'))
        eval(['clear ',a{i}]);
    end
end

goodUnits=xlsread('/mnt/muench_data/data/alexandra/MEA_data/sine_analysis/ListOfGoodFiltersPerND_sine.xls');
onOff_103=goodUnits(:,9);
goodUnits=goodUnits(:,1:8);

load(['/mnt/muench_data/data/alexandra/MEA_data/sine_analysis/HL10_accumParameters_NM12'])




% find cells with stable filters at each ND
figure
nds='8877665544332211';
celByND=zeros(16,103,4); % sequence: H30s, L30s, HL10 high, HL10 low
threshold=0.4;
for i=1:16
    subplot(4,4,i)
    hold on
    a=reshape(LinearFilter_H30s(:,1,i,:),500,103);
    plot(std(a(1:50,:))./std(a(85:200,:)),'.')
    celByND(i,1:103,1)=std(a(1:50,:))./std(a(85:200,:))<threshold;
    
    a=reshape(LinearFilter_L30s(:,1,i,:),500,103);
    plot(std(a(1:50,:))./std(a(85:200,:)),'.r')
    celByND(i,1:103,2)=std(a(1:50,:))./std(a(85:200,:))<threshold;
    
    a=reshape(mean(LinearFilter_HL10(:,1:2:end,i,:),2),500,103);
    plot(std(a(1:50,:))./std(a(85:200,:)),'.c')
    celByND(i,1:103,3)=std(a(1:50,:))./std(a(85:200,:))<threshold;
    
    a=reshape(mean(LinearFilter_HL10(:,2:2:end,i,:),2),500,103);
    plot(std(a(1:50,:))./std(a(85:200,:)),'.m')
    celByND(i,1:103,4)=std(a(1:50,:))./std(a(85:200,:))<threshold;
    
    axis([0 104 0 1.5])
    line([0 104],[threshold threshold],'color','k')
    title(['ND', nds(i),', red ON, blue OFF'])
    
end 


figure
nds='888777666555444333222111';
celByND_NM12=zeros(24,74,2); % sequence: HL10 high, HL10 low
cellNumber=74;
threshold=0.4;
for i=1:24
    subplot(4,6,i)
    hold on
    
    a=reshape(mean(LinearFilter(:,1:2:end,i,:),2),500,cellNumber);
    plot(std(a(1:50,:))./std(a(85:200,:)),'.b')
    celByND_NM12(i,1:cellNumber,1)=std(a(1:50,:))./std(a(85:200,:))<threshold;
    
    a=reshape(mean(LinearFilter(:,2:2:end,i,:),2),500,cellNumber);
    plot(std(a(1:50,:))./std(a(85:200,:)),'.r')
    celByND_NM12(i,1:cellNumber,2)=std(a(1:50,:))./std(a(85:200,:))<threshold;
    
    axis([0 cellNumber 0 1.5])
    line([0 cellNumber],[threshold threshold],'color','k')
    title(['ND', nds(i),', red ON, blue OFF'])
    
end 


%% Plots - for sine experiment only

%firing rate
figure
nds='87654321';
ranges=[-400 400 -400 400];
cnt=1;
for i=1:2:16
    subplot(2,4,cnt)
    hold on
    c=onOff_103<0;
    a=sum(mean(SpikeCount_HL10(i:i+1,c,1:2:end)),3)-sum(mean(SpikeCount_HL10(i:i+1,c,2:2:end)),3);
    b=mean(SpikeCount_H30s(i:i+1,c))-mean(SpikeCount_L30s(i:i+1,c));
    plot(a,b,'.b')
    [p,hoff(i)]=ttest(a-b);
    
    c=onOff_103>0;
    a=sum(mean(SpikeCount_HL10(i:i+1,c,1:2:end)),3)-sum(mean(SpikeCount_HL10(i:i+1,c,2:2:end)),3);
    b=mean(SpikeCount_H30s(i:i+1,c))-mean(SpikeCount_L30s(i:i+1,c));
    plot(a,b,'.r')
    [p,hon(i)]=ttest(a-b);    
    
    xlabel('d spike count HL10')
    ylabel('d spike count 30s')

    axis(ranges)
    line(ranges(1:2),[0 0],'color','k')
    line([0 0],ranges(3:4),'color','k')
    line(ranges(1:2),ranges(3:4),'color','k')
    title(['ND', nds(cnt),', red ON, blue OFF'])
    line([-50,-50],ranges(3:4),'color','k')
    line([50,50],ranges(3:4),'color','k')
    line(ranges(1:2),[-50,-50],'color','k')
    line(ranges(1:2),[50,50],'color','k')
    cnt=cnt+1;
end 



% firing rate at high vs firing rate delta
figure
nds='87654321';
ranges=[0 50 -400 400];
cnt=1;
for i=1:2:16
    subplot(2,4,cnt)
    hold on
    c=onOff_103<0;
    a=mean(SpikeCount_H30s(i:i+1,c))/30; % roughly Hz
    b=mean(SpikeCount_H30s(i:i+1,c))-mean(SpikeCount_L30s(i:i+1,c));
    plot(a,b,'.b')
    [p,hoff(i)]=ttest(a-b);
    
    c=onOff_103>0;
    a=mean(SpikeCount_H30s(i:i+1,c))/30; % roughly Hz
    b=mean(SpikeCount_H30s(i:i+1,c))-mean(SpikeCount_L30s(i:i+1,c));
    plot(a,b,'.r')
    [p,hon(i)]=ttest(a-b);    
    
    xlabel('base firing rate at HL10 high')
    ylabel('d firing rate 30s')

    axis(ranges)
    line(ranges(1:2),[0 0],'color','k')
    line([0 0],ranges(3:4),'color','k')
    title(['ND', nds(cnt),', red ON, blue OFF'])
%     line([-50,-50],ranges(3:4),'color','k')
%     line([50,50],ranges(3:4),'color','k')
%     line(ranges(1:2),[-50,-50],'color','k')
%     line(ranges(1:2),[50,50],'color','k')
    cnt=cnt+1;
end 


% firing rate at high vs firing rate at low
figure
nds='87654321';
ranges=[0 60 0 60];
cnt=1;
for i=1:2:16
    subplot(2,4,cnt)
    hold on
    c=goodUnits(:,cnt)>0&onOff_103<0;
    a=mean(SpikeCount_H30s(i:i+1,c))/30; % roughly Hz
    b=mean(SpikeCount_L30s(i:i+1,c))/30;
    plot(a,b,'.b')
    [p,hoff(i)]=ttest(a-b);
    
    c=goodUnits(:,cnt)>0&onOff_103>0;
    a=mean(SpikeCount_H30s(i:i+1,c))/30; % roughly Hz
    b=mean(SpikeCount_L30s(i:i+1,c))/30;
    plot(a,b,'.r')
    [p,hon(i)]=ttest(a-b);    
    
    xlabel('spike count at high')
    ylabel('spike count at low')

    axis(ranges)
    line(ranges(1:2),[0 0],'color','k')
    line([0 0],ranges(3:4),'color','k')
    title(['ND', nds(cnt),', red ON, blue OFF'])
    line(ranges(1:2),ranges(3:4),'color','k')
    cnt=cnt+1;
end 


%zero crossing
figure
nds='87654321';
ranges=[-40 40 -40 40];
cnt=1;
for i=1:2:16
    subplot(2,4,cnt)
    hold on
    c=goodUnits(:,cnt)==1&onOff<0;
    a=mean(ind_HL10(1:2:end,i,c),1)-mean(ind_HL10(2:2:end,i,c),1);
    a1=mean(ind_HL10(1:2:end,i+1,c),1)-mean(ind_HL10(2:2:end,i+1,c),1);
    a=mean([a; a1]);
    a=reshape(a,size(a,3),1);
    
    b=ind_H30s(1,i,c)-ind_L30s(1,i,c);
    b1=ind_H30s(1,i+1,c)-ind_L30s(1,i+1,c);
    b=mean([b; b1]);
    b=reshape(b,size(b,3),1);
    plot(a,b,'.b')    
    [~,hoff(i)]=ttest(a-b);

    c=goodUnits(:,cnt)==1&onOff>0;
    a=mean(ind_HL10(1:2:end,i,c),1)-mean(ind_HL10(2:2:end,i,c),1);
    a1=mean(ind_HL10(1:2:end,i+1,c),1)-mean(ind_HL10(2:2:end,i+1,c),1);
    a=mean([a; a1]);
    a=reshape(a,size(a,3),1);
    
    b=ind_H30s(1,i,c)-ind_L30s(1,i,c);
    b1=ind_H30s(1,i+1,c)-ind_L30s(1,i+1,c);
    b=mean([b; b1]);
    b=reshape(b,size(b,3),1);
    plot(a,b,'.b') 
    plot(a,b,'.r')
    [~,hoff(i)]=ttest(a-b);
    
    
    xlabel('zero crossing HL10')
    ylabel('zero crossing 30s')
    axis(ranges)
    line(ranges(1:2),[0 0],'color','k')
    line([0 0],ranges(3:4),'color','k')
    line(ranges(1:2),ranges(3:4),'color','k')
    title(['ND', nds(cnt),', red ON, blue OFF'])
    cnt=cnt+1;
end 

% zero crossing at high vs zero crossing  at low
figure
nds='87654321';
ranges=[0 300 0 300];
cnt=1;
for i=1:2:16
    subplot(2,4,cnt)
    hold on
    c=goodUnits(:,cnt)>0&onOff_103<0;
    a=mean(ind_H30s(1:2:end,i:i+1,c));
    b=mean(ind_L30s(1:2:end,i:i+1,c));
    a=reshape(a,size(a,3),1);
    b=reshape(b,size(b,3),1);
    plot(a,b,'.b')
    [p,hoff(cnt)]=ttest(a-b);
    
    c=goodUnits(:,cnt)>0&onOff_103>0;
    a=mean(ind_H30s(1:2:end,i:i+1,c));
    b=mean(ind_L30s(1:2:end,i:i+1,c));
    a=reshape(a,size(a,3),1);
    b=reshape(b,size(b,3),1);
    plot(a,b,'.r')
    [p,hon(cnt)]=ttest(a-b);    
    
    xlabel('zero x at high')
    ylabel('zero x at low')

    axis(ranges)
    line(ranges(1:2),[0 0],'color','k')
    line([0 0],ranges(3:4),'color','k')
    title(['ND', nds(cnt),', red ON, blue OFF'])
    line(ranges(1:2),ranges(3:4),'color','k')
    cnt=cnt+1;
end 



%firing rate vs zero crossing
figure
nds='87654321';
ranges=[-700 700 -50 50]*2;
cnt=1;
for i=1:2:16
    subplot(2,4,cnt)
    hold on
    c=goodUnits(:,cnt)>0&onOff_103<0;
    a=mean(SpikeCount_H30s(i:i+1,c))-mean(SpikeCount_L30s(i:i+1,c));    
    b=mean(ind_H30s(1,i:i+1,c))-mean(ind_L30s(1,i:i+1,c));
    b=reshape(b,1,size(b,3));
    plot(a,b,'.b')    

    c=goodUnits(:,cnt)>0&onOff_103>0;
    a=mean(SpikeCount_H30s(i:i+1,c))-mean(SpikeCount_L30s(i:i+1,c));    
    b=mean(ind_H30s(1,i:i+1,c))-mean(ind_L30s(1,i:i+1,c));
    b=reshape(b,1,size(b,3));
    plot(a,b,'.r')
    
    
    xlabel('spike count H30s-L30s')
    ylabel('zero crossing H30s-L30s')
    axis(ranges)
    line(ranges(1:2),[0 0],'color','k')
    line([0 0],ranges(3:4),'color','k')
    title(['ND', nds(cnt),', red ON, blue OFF'])
    cnt=cnt+1;
end 


%peak
figure
nds='8877665544332211';
ranges=[-200 200 -200 200];
for i=1:16
    subplot(4,4,i)
    hold on
    c=celByND(i,:,1)&celByND(i,:,2)&celByND(i,:,3)&celByND(i,:,4)&onOff_H30s<0;
    a=sum(peak_HL10(1:2:end,i,c),1)-sum(peak_HL10(2:2:end,i,c),1);
    a=reshape(a,size(a,3),1);
    
    b=peak_H30s(1,i,c)-peak_L30s(1,i,c);
    b=reshape(b,size(b,3),1);
    plot(a,b,'.b')    
    [~,hoff(i)]=ttest(a-b);

    c=celByND(i,:,1)&celByND(i,:,2)&celByND(i,:,3)&celByND(i,:,4)&onOff_H30s>0;
    a=sum(peak_HL10(1:2:end,i,c),1)-sum(peak_HL10(2:2:end,i,c),1);
    a=reshape(a,size(a,3),1);
    
    b=peak_H30s(1,i,c)-peak_L30s(1,i,c);
    b=reshape(b,size(b,3),1);
    plot(a,b,'.r')
    [~,hoff(i)]=ttest(a-b);
    
    
    xlabel('HL10')
    ylabel('30s')
    axis(ranges)
    line(ranges(1:2),[0 0],'color','k')
    line([0 0],ranges(3:4),'color','k')
    line(ranges(1:2),ranges(3:4),'color','k')
    title(['ND', nds(i),', red ON, blue OFF'])
end 


%peak latency
figure
nds='8877665544332211';
ranges=[-20 20 -20 20];
for i=1:16
    subplot(4,4,i)
    hold on
    c=celByND(i,:,1)&celByND(i,:,2)&celByND(i,:,3)&celByND(i,:,4)&onOff_H30s<0;
    a=mean(peakLatency_HL10(1:2:end,i,c),1)-mean(peakLatency_HL10(2:2:end,i,c),1);
    a=reshape(a,size(a,3),1);
    
    b=peakLatency_H30s(1,i,c)-peakLatency_L30s(1,i,c);
    b=reshape(b,size(b,3),1);
    plot(a,b,'.b')    
    [~,hoff(i)]=ttest(a-b);

    c=celByND(i,:,1)&celByND(i,:,2)&celByND(i,:,3)&celByND(i,:,4)&onOff_H30s>0;
    a=mean(peakLatency_HL10(1:2:end,i,c),1)-mean(peakLatency_HL10(2:2:end,i,c),1);
    a=reshape(a,size(a,3),1);
    
    b=peakLatency_H30s(1,i,c)-peakLatency_L30s(1,i,c);
    b=reshape(b,size(b,3),1);
    plot(a,b,'.r')
    [~,hoff(i)]=ttest(a-b);
    
    
    xlabel('HL10')
    ylabel('30s')
    axis(ranges)
    line(ranges(1:2),[0 0],'color','k')
    line([0 0],ranges(3:4),'color','k')
    line(ranges(1:2),ranges(3:4),'color','k')
    title(['ND', nds(i),', red ON, blue OFF'])
end 


%gain
figure
nds='8877665544332211';
ranges=[0 0.05 0 0.05];
for i=1:16
    subplot(4,4,i)
    hold on
    c=celByND(i,:,1)&celByND(i,:,2)&celByND(i,:,3)&celByND(i,:,4)&onOff_H30s<0;
    a=mean(nonLinear_HL10(2,1:2:end,i,c),2)./mean(nonLinear_HL10(2,2:2:end,i,c),2);
    a=reshape(a,size(a,4),1);
    
    b=nonLinear_H30s(2,1,i,c)./nonLinear_L30s(2,1,i,c);
    b=reshape(b,size(b,4),1);
    plot(a,b,'.b')    
    [~,hoff(i)]=ttest(a-b);

    c=celByND(i,:,1)&celByND(i,:,2)&celByND(i,:,3)&celByND(i,:,4)&onOff_H30s>0;
    a=mean(nonLinear_HL10(2,1:2:end,i,c),2)./mean(nonLinear_HL10(2,2:2:end,i,c),2);
    a=reshape(a,size(a,4),1);
    
    b=nonLinear_H30s(2,1,i,c)./nonLinear_L30s(2,1,i,c);
    b=reshape(b,size(b,4),1);
    plot(a,b,'.r')
    [~,hoff(i)]=ttest(a-b);
    
    
    xlabel('HL10')
    ylabel('30s')
    axis(ranges)
    line(ranges(1:2),[0 0],'color','k')
    line([0 0],ranges(3:4),'color','k')
    line(ranges(1:2),ranges(3:4),'color','k')
    title(['ND', nds(i),', red ON, blue OFF'])
end 


%adjusted gain
figure
nds='8877665544332211';
ranges=[0 0.01 0 0.01];
for i=1:16
    subplot(4,4,i)
    hold on
    c=celByND(i,:,1)&celByND(i,:,2)&celByND(i,:,3)&celByND(i,:,4)&onOff_H30s<0;
    a=mean(nonLinear_adjusted_HL10(1:2:end,i,c),1)./mean(nonLinear_adjusted_HL10(2:2:end,i,c),1);
    a=reshape(a,size(a,3),1);
    
    b=nonLinear_adjusted_H30s(1,i,c)./nonLinear_adjusted_L30s(1,i,c);
    b=reshape(b,size(b,3),1);
    plot(a,b,'.b')    
    [~,hoff(i)]=ttest(a-b);

    c=celByND(i,:,1)&celByND(i,:,2)&celByND(i,:,3)&celByND(i,:,4)&onOff_H30s>0;
    a=mean(nonLinear_adjusted_HL10(1:2:end,i,c),1)./mean(nonLinear_adjusted_HL10(2:2:end,i,c),1);
    a=reshape(a,size(a,3),1);
    
    b=nonLinear_adjusted_H30s(1,i,c)./nonLinear_adjusted_L30s(1,i,c);
    b=reshape(b,size(b,3),1);
    plot(a,b,'.r')
    [~,hoff(i)]=ttest(a-b);
    
    
    xlabel('HL10')
    ylabel('30s')
    axis(ranges)
    line(ranges(1:2),[0 0],'color','k')
    line([0 0],ranges(3:4),'color','k')
    line(ranges(1:2),ranges(3:4),'color','k')
    title(['ND', nds(i),', red ON, blue OFF'])
end 


%% Plots - sine and NM12


%firing rate vs zero crossing
figure
nds='87654321';
ranges=[-250 250 -150 150];
cnt_sine=1;
cnt_nm=1;
for i=1:8
    subplot(3,3,i)
    hold on
    % sine
    tmp=logical(sum(celByND(cnt_sine:cnt_sine+1,:,3)));
    tmp1=logical(sum(celByND(cnt_sine:cnt_sine+1,:,4)));
    c=tmp&tmp1&onOff_HL10<0;
    a=mean(SpikeCount_HL10(cnt_sine,c,1:2:end),3)-mean(SpikeCount_HL10(cnt_sine,c,2:2:end),3); 
    a1=mean(SpikeCount_HL10(cnt_sine+1,c,1:2:end),3)-mean(SpikeCount_HL10(cnt_sine+1,c,2:2:end),3);   
    a=mean([a;a1]);
    
    a=mean(SpikeCount_HL10(cnt_sine,c,1:2:end),3)-mean(SpikeCount_HL10(cnt_sine,c,2:2:end),3); 
    a1=mean(SpikeCount_HL10(cnt_sine+1,c,1:2:end),3)-mean(SpikeCount_HL10(cnt_sine+1,c,2:2:end),3);   
    a=mean([a;a1]);
    
    b=mean(ind_HL10(1:2:end,cnt_sine,c))-mean(ind_HL10(2:2:end,cnt_sine,c));
    b1=mean(ind_HL10(1:2:end,cnt_sine+1,c))-mean(ind_HL10(2:2:end,cnt_sine+1,c));
    b=reshape(b,1,size(b,3));
    b1=reshape(b1,1,size(b1,3));
    b=mean([b;b1]);
    
    plot(a,b,'.b')    


    tmp=logical(sum(celByND(cnt_sine:cnt_sine+1,:,3)));
    tmp1=logical(sum(celByND(cnt_sine:cnt_sine+1,:,4)));
    c=tmp&tmp1&onOff_HL10>0;
    a=mean(SpikeCount_HL10(cnt_sine,c,1:2:end),3)-mean(SpikeCount_HL10(cnt_sine,c,2:2:end),3); 
    a1=mean(SpikeCount_HL10(cnt_sine+1,c,1:2:end),3)-mean(SpikeCount_HL10(cnt_sine+1,c,2:2:end),3);   
    a=mean([a;a1]);
    
    a=mean(SpikeCount_HL10(cnt_sine,c,1:2:end),3)-mean(SpikeCount_HL10(cnt_sine,c,2:2:end),3); 
    a1=mean(SpikeCount_HL10(cnt_sine+1,c,1:2:end),3)-mean(SpikeCount_HL10(cnt_sine+1,c,2:2:end),3);   
    a=mean([a;a1]);
    
    b=mean(ind_HL10(1:2:end,cnt_sine,c))-mean(ind_HL10(2:2:end,cnt_sine,c));
    b1=mean(ind_HL10(1:2:end,cnt_sine+1,c))-mean(ind_HL10(2:2:end,cnt_sine+1,c));
    b=reshape(b,1,size(b,3));
    b1=reshape(b1,1,size(b1,3));
    b=mean([b;b1]);
    
    plot(a,b,'.r')
    
    cnt_sine=cnt_sine+2;
    
    % nm12
    tmp=logical(sum(celByND_NM12(cnt_nm:cnt_nm+2,:,1)));
    tmp1=logical(sum(celByND_NM12(cnt_nm:cnt_nm+2,:,2)));
    c=tmp&tmp1&onOff<0;
    a=mean(SpikeCount(cnt_nm,c,1:2:end),3)-mean(SpikeCount(cnt_nm,c,2:2:end),3); 
    a1=mean(SpikeCount(cnt_nm+1,c,1:2:end),3)-mean(SpikeCount(cnt_nm+1,c,2:2:end),3);   
    a2=mean(SpikeCount(cnt_nm+2,c,1:2:end),3)-mean(SpikeCount(cnt_nm+2,c,2:2:end),3);
    a=mean([a;a1;a2]);
    
    
    b=mean(ind(1:2:end,cnt_nm,c))-mean(ind(2:2:end,cnt_nm,c));
    b1=mean(ind(1:2:end,cnt_nm+1,c))-mean(ind(2:2:end,cnt_nm+1,c));
    b2=mean(ind(1:2:end,cnt_nm+2,c))-mean(ind(2:2:end,cnt_nm+2,c));
    b=mean([b;b1;b2],1);
    b=reshape(b,1,size(b,3));    
    
    plot(a,b,'.c')    

    tmp=logical(sum(celByND_NM12(cnt_nm:cnt_nm+2,:,1)));
    tmp1=logical(sum(celByND_NM12(cnt_nm:cnt_nm+2,:,2)));
    c=tmp&tmp1&onOff>0;
    a=mean(SpikeCount(cnt_nm,c,1:2:end),3)-mean(SpikeCount(cnt_nm,c,2:2:end),3); 
    a1=mean(SpikeCount(cnt_nm+1,c,1:2:end),3)-mean(SpikeCount(cnt_nm+1,c,2:2:end),3);   
    a2=mean(SpikeCount(cnt_nm+2,c,1:2:end),3)-mean(SpikeCount(cnt_nm+2,c,2:2:end),3);
    a=mean([a;a1;a2]);
    
    
    b=mean(ind(1:2:end,cnt_nm,c))-mean(ind(2:2:end,cnt_nm,c));
    b1=mean(ind(1:2:end,cnt_nm+1,c))-mean(ind(2:2:end,cnt_nm+1,c));
    b2=mean(ind(1:2:end,cnt_nm+2,c))-mean(ind(2:2:end,cnt_nm+2,c));
    b=mean([b;b1;b2],1);
    b=reshape(b,1,size(b,3)); 
    
    plot(a,b,'.m')

    cnt_nm=cnt_nm+3;
    
    
    xlabel('firing rate HL10 high - HL10 low')
    ylabel('zero crossing HL10 high - HL10 low')
    axis(ranges)
    line(ranges(1:2),[0 0],'color','k')
    line([0 0],ranges(3:4),'color','k')
    title(['ND', nds(i),', red ON, blue OFF'])
end 

figure
for i=1:24
    subplot(4,6,i)
    tmp=reshape(SpikeCount(i,:,:),74,6)';
    for j=1:74
        x(i,j)=corr(tmp(:,j),sigALL');
    end
    plot(x(i,onOff>0),'*r')
    hold on
    plot(x(i,onOff<0),'*b')
end