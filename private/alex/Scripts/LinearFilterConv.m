
clear
date='20121023';
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date])
load([path2save,'protocols_',codeWord])

filter_length=500;

if length(unique(startTimes))>1
    for i=1:size(correctedProtocols,3)
        correctedProtocols(:,1,i)=correctedProtocols(:,1,i)-startTimes(i);
    end
else
    correctedProtocols(:,1,:)=correctedProtocols(:,1,:)-startTimes(1);
end
LinearFilterConv=zeros(700,192,39);
for cnt=1:length(units)
    
    load([mainpath,'units/',units(cnt).name]);
    for i=1:length(file_list)
        spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
        if ~isempty(spikes)
            spikes=spikes-startTimes(i)+1;
            for t=1:dim
                startPoint=correctedProtocols(per*(t-1)+1,1,i)+1;
                endPoint=correctedProtocols(per*(t-1)+dur,1,i);
                spikesTMP=spikes(spikes>startPoint+500&spikes<endPoint);
                convSpikes=convolved(spikesTMP,40,size(flicker,1));
                convSpikes=convSpikes(121:end-120);
%                 convSpikes=convSpikes-mean(convSpikes(200:1700));
                [a,b]=hist(convSpikes,20);
                clear lfc cfc a
                cfc=zeros(20,1);
                i=49
                tmp1=0;
                tmp2=0;
                for kk=1:length(b)-1
                    tmp=0;                    
                    for j=201:length(flicker)-1000
                        if convSpikes(j+499)>=b(kk)&&convSpikes(j+499)<b(kk+1)
                            tmp(1:700)=tmp+flicker(j:j+699,i)'*convSpikes(j+499);
                            tmp1(1:700)=tmp1+flicker(j:j+699,i)'*convSpikes(j+499);
                            if kk>7
                                tmp2(1:700)=tmp2+flicker(j:j+699,i)'*convSpikes(j+499);
                            end
                            cfc(kk)=cfc(kk)+1;
                        end
                    end
                    lfc(700:-1:1,kk)=tmp/(sum(abs(tmp)));
                end
                a(700:-1:1)=tmp1/(sum(abs(tmp1)));
                b(700:-1:1)=tmp2/(sum(abs(tmp2)));
                figure
                for i=1:19
                    subplot(5,4,i)
                    plot(lfc(:,i))
                    hold on
                    plot(a,'r')
                    plot(b,'g')
                    title([int2str(cfc(i)), '   ', num2str(b(i))])
                end
                figure;plot(lfc(:,1:5))
                figure;plot(lfc(:,5:10))
                
                    
%                 n=zeros(length(spikesTMP),filter_length);
%                 for k=1:filter_length
%                     n(:,k)=flicker(spikesTMP-k+1,i);
%                 end
%                 figure
%                plot(n')
            end
        end

    end
end



seq1=zeros(55000,39);
for cnt=1:39
    
    load([mainpath,'units/',units(cnt).name]);
    for i=49:6:72
        spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
        spikes=spikes-startTimes(i)+1;
        startPoint=correctedProtocols(1,1,i)+1;
        endPoint=correctedProtocols(dur,1,i);
        spikesTMP=spikes(spikes>startPoint+500&spikes<endPoint);
        convSpikes=convolved(spikesTMP,40,size(flicker,1));
        convSpikes=convSpikes(121:end-120);
        seq1(1:55000,cnt)=seq1(1:55000,cnt)+convSpikes(1:55000)';
    end
end

seq1=seq1/3;
seq1(:,1)=seq1(:,1)/max(seq1(:,1));
seq1(:,2)=seq1(:,2)/max(seq1(:,2));
seq1(:,3)=seq1(:,3)/max(seq1(:,3));
figure
plot(seq1)

clear a
c=1;
for i=1:10:54699
    a(c)=corr(seq1(i:i+200,1),seq1(i:i+200,2));
    c=c+1;
end

tst=1:10:54699;
b=find(a<-0.5);
tr=[];
for k=1:length(b)-1
    if b(k+1)-b(k)>10
        tr=[tr b(k)];
    end
end
figure
for i=1:length(tr)
    subplot(8,8,i)
    plot(seq1(tst(tr(i)):tst(tr(i))+200,1))
    hold on
    plot(seq1(tst(tr(i)):tst(tr(i))+200,2),'r')
end

[ind lags]=xcorr(seq1(:,1),seq1(:,2),'coeff');



figure;imagesc(corr(seq1))

figure;imagesc(corr(seq1(:,onOff<0)))

figure;imagesc(corr(seq1(:,onOff>0)))

k=find(onOff>0);
k(wrongCells)

figure;imagesc(corr(seq1(:,k(~wrongCells))))
