
datarun=datarunD;
datarun=load_sta(datarun,'load_sta','all');
clear dtarunD

m=[931,1471];%,1276]
sta=zeros(320,320,3);
for i=1:length(m)
    cellID=find(datarun.cell_ids==m(i));
    
    a=datarun.stas.stas{cellID, 1};
    a=squeeze(a);
    a=a(:,:,5);    
    sta(:,:,i)=a;
end
figure
imagesc(sta*5+0.5)

k=reshape(sta(:,:,1:2),320*320,2);

m=datarun.cell_types{1, 3}.cell_ids(1:2)';
sta=zeros(320,320);
for i=1:length(m)
    cellID=find(datarun.cell_ids==m(i));
    
    a=datarun.stas.stas{cellID, 1};
    a=squeeze(a);
    a=a(:,:,5);    
    sta(:,:)=sta+a;
end
figure
colormap gray
imagesc(sta)

m=[151,166];%,1276]
sta=zeros(320,320,3);
for i=1:length(m)
    cellID=find(datarun.cell_ids==m(i));
    
    a=datarun.stas.stas{cellID, 1};
    a=squeeze(a);
    a=a(:,:,5);    
    sta(:,:,i)=a;
end
figure
imagesc(sta*3+0.5)



a=reshape(a,320*320,6);

c=corr(a');

figure
colormap gray
imagesc(a(:,:,5))

fd=find_layerers(datarun, 64)

for i=1:length(fd)
    datarun.cell_ids(fd(i))
end

temp_cell_ID = m(1);
temp_rf = get_rf(datarun, temp_cell_ID);
sig_stix = significant_stixels(temp_rf, 'thresh', 5);
[timecourse, params] = time_course_from_sta(sta, sig_stix);
figure
plot(timecourse)

a=datarun.stas.stas{temp_cell_ID, 1};
a=squeeze(a);
a=a(:,:,2:end);
figure
plot(reshape(a(120,146,:),5,1))
imagesc(a(:,:,4))


t=reshape(a,320*320,5);
r=repmat(timecourse,320*320,1);
k=t.*r;

k=reshape(k,320,320,5);

figure;
colormap gray
imagesc(k(:,:,4))


t=reshape(a,320*320,5);
b=zeros(length(t),1);
for i=1:length(t)
    b(i)=corr(timecourse',t(i,:)');
end

c=reshape(b,320,320);
figure;
colormap gray
imagesc(c)
figure;
colormap gray
imagesc(a(:,:,4))

k=a(:,:,4).*c;
figure;
colormap gray
imagesc(k)

robust_std(k(:))
b=a(:,:,4);
robust_std(b(:))
robust_std(c(:))

max(k(:))/robust_std(k(:))

max(b(:))/robust_std(b(:))





cellID=find(datarun.cell_ids==151);
spikes=datarun.spikes{cellID};

cellID1=find(datarun.cell_ids==166);
spikes1=datarun.spikes{cellID1};



% spikes are in s. convert to ms
spikes=round(spikes*1000);
spikes(spikes<5*refresh*2)=[];
spikes1=round(spikes1*1000);
spikes1(spikes1<5*refresh*2)=[];
%fr - frames - are in ms
fr=round(triggers(1)*1000:refresh:triggers(end)*1000);
timethresh=3;
sta=zeros(height,width,5); %height, width, frames back
tic
nex=0.1;
cnt=1;
for i=spikes'
    if min(abs(spikes1-i))<timethresh
        cnt=cnt+1;
        start=find(fr>i,1)-7;
        for j=1:5
            F = round(mvi.getFrame(start+j).getBuffer);
            sta(:,:,j) = sta(:,:,j) + reshape(F(1:3:end),width,height)';
        end
        if i>spikes(end)*nex
            disp(nex)
            nex=nex+0.1;
        end
    end
end
toc
figure
colormap gray
imagesc(sta(:,:,4))

a=sta(:,:,4);









% get indices for cells of interest.
cellspec = {1,2,3,4};
temp_indices = get_cell_indices(datarun, cellspec);
numRGCs = length(temp_indices);
 timecourse=zeros(6,temp_indices(end));
all_sig_stix = zeros(datarun.stimulus.field_height*datarun.stimulus.field_width, temp_indices(end));
% process one cell at a time
for cc = 1:numRGCs
    % get cell ID and RF
    temp_cell_ID = datarun.cell_ids(temp_indices(cc));
    temp_rf = get_rf(datarun, temp_cell_ID);        
    
    % check to see if stim is RGB
    sig_stix = significant_stixels(temp_rf, 'thresh', 5);
    
    if ~isempty(temp_rf)&&nnz(sig_stix)
        
        all_sig_stix(:,temp_indices(cc)) = sig_stix(:);
        sta=datarun.stas.stas{temp_indices(cc), 1};
        [timecourse_tmp,~] = time_course_from_sta(sta, sig_stix);
        timecourse(:,temp_indices(cc))=timecourse_tmp;
    end
end

cell_list=cell(1,temp_indices(end));

for j=1:numRGCs
    cellID=temp_indices(j);
    tmp=find(all_sig_stix(:,cellID));
    for i=1:length(tmp)
        cell_list{cellID}=[cell_list{cellID} find(all_sig_stix(tmp(i),:))];
    end
    cell_list{cellID}=unique(cell_list{cellID});
    cell_list{cellID}(cell_list{cellID}==cellID)=[];   
end


for j=1:numRGCs
    cellID=temp_indices(j);
    a=round(datarun.stas.fits{cellID}.mean);
    y=[max(1,a(1)-50) min(320,a(1)+50)];
    x=[max(1,a(2)-50) min(320,a(2)+50)];
    sta=datarun.stas.stas{cellID, 1};
    sta=squeeze(sta);
%     sta=sta(x(1):x(2),y(1):y(2),:);
    dim=size(sta);
    
    figure
    subplot(1,3,1)
    colormap gray
    imagesc(sta(:,:,5))
    
    sta1=reshape(sta,prod(dim(1:2)),6);
    t=timecourse(:,cellID);
    a=corr(sta1',t);
    a=reshape(a,dim(1),dim(2));

    subplot(1,3,2)
    colormap gray
    imagesc(a)
    
    subplot(1,3,3)
    k=sta(:,:,5).*a;
    colormap gray
    imagesc(k)
end

 sig_stix = significant_stixels(k, 'thresh', 12);
 sig_stix=full(sig_stix);
 figure
 subplot(1,2,1)
 imagesc(sig_stix)
 sig_stix = significant_stixels(sta(:,:,5), 'thresh', 5);
 sig_stix=full(sig_stix);
 subplot(1,2,2)
 imagesc(sig_stix)
 
 
 
 
 
 stixel_threshold = gen_params.stixel_threshold; 

% get indices for cells of interest.
temp_indices = get_cell_indices(datarun, cellspec);
numRGCs = length(temp_indices);
 
% initialze a big matrix
normalized_RFs = zeros(numRGCs, prod([size(datarun.stas.rfs{temp_indices(1)},1), size(datarun.stas.rfs{temp_indices(1)},2)]));
all_sig_stix = zeros(prod([size(datarun.stas.rfs{temp_indices(1)},1), size(datarun.stas.rfs{temp_indices(1)},2)]), numRGCs);
% process one cell at a time
for cc = 1:numRGCs
    % get cell ID and RF
    temp_cell_ID = datarun.cell_ids(temp_indices(cc));
%     temp_rf = get_rf(datarun, temp_cell_ID);        
    
    % check to see if stim is RGB
    if size(temp_rf, 3) == 3
       
        % check to see if SBC
        if ~isempty(find(datarun.cell_types{5}.cell_ids == temp_cell_ID))
            
            % get just the blue gun
            temp_rf = squeeze(temp_rf(:,:,3));
            
            % get sig stixel
            sig_stix = significant_stixels(temp_rf, 'thresh', 5);
            
%             figure(1)
%             imagesc(temp_rf); colormap gray; 
%             pause
%             imagesc(sig_stix)
%             pause
        else
            % otherwise sum all guns
            temp_rf = sum(temp_rf, 3);

            sig_stix = significant_stixels(temp_rf, 'thresh', stixel_threshold);
        end
    else
        % if already BW, calculate significant stixels
                
        cellID=temp_indices(cc);
        sta=datarun.stas.stas{cellID, 1};
        sta=squeeze(sta);
        dim=size(sta);
        
        sta1=reshape(sta,prod(dim(1:2)),6);
        t=timecourse(:,cellID);
        a=corr(sta1',t);
        a=reshape(a,dim(1),dim(2));
        k=sta(:,:,5).*a;
        temp_rf=k;
    
        sig_stix = significant_stixels(temp_rf, 'thresh', 15);
    end
    
    if ~isempty(temp_rf)
        % normalize by the SNR
        temp_N = robust_std(reshape(temp_rf, 1,[]), 5);     % estimate the noise
        temp_S = mean(temp_rf(sig_stix));     % estimate the signal
        temp_SNR = temp_S/temp_N^2; % estimate the SNR
        
        % threshold the RF and scale by SNR
        thresh_rf = zeros(size(temp_rf));
        thresh_rf(sig_stix) = temp_rf(sig_stix) * temp_SNR;
        normalized_RFs(cc,:) = reshape(thresh_rf, 1,[]);
        
        % Save for refinement stage
        
        all_sig_stix(:,cc) = sig_stix(:);
    end
end



%%

clear
datarun.names.rrs_neurons_path='/Volumes/Analysis/2010-03-05-2/data001/data001.neurons';
    
opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',true);
datarun=load_data(datarun,opt);

triggers=datarun.triggers; %onsets of the stimulus presentation

mdf_file='/Volumes/Analysis/deprecated/movie-xml2/BW-1-8-0.48-11111-320x320.xml';

[mov,height,width,duration,refresh] = get_movie_ath('/Volumes/Analysis/deprecated/movie-xml2/BW-1-8-0.48-11111-320x320.xml',...
    triggers,2,6);

[mvi] = load_movie(mdf_file, triggers);

figure
colormap gray
imagesc(mov(:,:,2))

cellID=find(datarun.cell_ids==151);


spikes=datarun.spikes{cellID};


cellID=find(datarun.cell_ids==1471);

spikes=sort([spikes; datarun.spikes{cellID}]);


% spikes are in s. convert to ms
spikes=round(spikes*1000);
spikes(spikes<5*refresh*2)=[];
%fr - frames - are in ms
fr=round(triggers(1)*1000:refresh:triggers(end)*1000);

sta=zeros(height,width,5); %height, width, frames back
tic
nex=0.1;
for i=spikes'
 
    start=find(fr>i,1)-7; 
    for j=1:5
        F = round(mvi.getFrame(start+j).getBuffer);
        sta(:,:,j) = sta(:,:,j) + reshape(F(1:3:end),width,height)';
    end
    if i>spikes(end)*nex
        disp(nex)
        nex=nex+0.1;
    end
end
toc
figure
colormap gray
imagesc(sta(:,:,4))






spikes=datarun.spikes{cellID};

% spikes are in s. convert to ms
spikes=round(spikes*1000);
spikes(spikes<5*refresh*2)=[];
spikes(end-200:end)=[];
%fr - frames - are in ms
fr=round(triggers(1)*1000:refresh:triggers(end)*1000);

sta=zeros(height,width,5); %height, width, frames back
tic
nex=0.1;
spikes=spikes(2:2:end);
for i=spikes'
 
    start=find(fr>i,1)-7; 
    for j=1:5
        F = round(mvi.getFrame(start+j).getBuffer);
        sta(:,:,j) = sta(:,:,j) + reshape(F(1:3:end),width,height)';
    end
    if i>spikes(end)*nex
        disp(nex)
        nex=nex+0.1;
    end
end
toc
figure
colormap gray
imagesc(sta(:,:,4))

axis([80 115 100 145])

sta1=sta;

sta1=sta1(80:115,100:145,:);
sta=sta(80:115,100:145,:);


a=reshape(sta1,36*46,5);
b=reshape(sta,36*46,5);
for i=1:1656
        t(i)=corr(a(i,2:4)',b(i,2:4)');
end

t=reshape(t,36,46);
figure
colormap gray
imagesc(t)


sta=zeros(100,100,4); %height, width, frames back
codeWords=zeros(10000,14000);
tmp=zeros(10000,4);
tic
nex=0.01;
cnt=1;
for i=spikes'
    if cnt>14000
        break;
    end
        start=find(fr>i,1)-6;
        start

        for j=1:4
            F = round(mvi.getFrame(start+j).getBuffer);
            F=reshape(F(1:3:end),width,height)';
            F=F(51:150,51:150);
            sta(:,:,j) = sta(:,:,j) + F;
            tmp(:,j)=F(:);
        end
        tmp=bin2dec(int2str(tmp));
        codeWords(:,cnt)=tmp;
        %     for j=0:15
        %         codeWords(a==j,j+1)=codeWords(a==j,j+1)+1;
        %     end
            cnt

        cnt=cnt+1;
    
%     sta = sta + reshape(tmp,40,40,4);
%     tmp2=tmp2+tmp;
end
toc

codeWords=codeWords(:,1:cnt);

a=codeWords(1,:);
b=zeros(1,16);
for i=0:15
    b(i+1)=sum(a==i);
end

figure
bar(b)

figure
for i=1:16
    a=codeWords;
    a(codeWords~=i-1)=0;
    a(codeWords==i-1)=1;
    
    a=sum(a,2);
    subplot(4,4,i)
    colormap gray
    imagesc(reshape(a,100,100));
end


a=reshape(F,50,50)';
figure;
imagesc(a)


save('/Users/alexth/Desktop/codeWords4_60.mat','codeWords')


codeWords(1:10,:)
figure
for i=1:16
    subplot(4,4,i)
    colormap('gray')
    imagesc(reshape(codeWords(:,i),50,50))
    title(int2str(mean(codeWords(:,i))))
end

figure
clear a b
for i=0:15
   subplot(4,4,i+1)
   colormap('gray')
   a=dec2bin(i,4);
   for j=1:4
       b(j)=eval(a(j));
   end
   imagesc(b)
   axis([0.5 4.5 0.5 1.5])
   title([int2str(mean(codeWords(:,i+1))),'+-', int2str(std(codeWords(:,i+1)))])
end



a=codeWords(:,2)+codeWords(:,3)+codeWords(:,8)+codeWords(:,11)+codeWords(:,12);
b=codeWords(:,5)+codeWords(:,6)+codeWords(:,9)+codeWords(:,13)+codeWords(:,14);
figure
colormap('gray')
imagesc(reshape(a-b,50,50))

figure
colormap('gray')
c=1-std(codeWords,0,2);
imagesc(reshape(c,60,60))
d=sub2ind([60,60],40,26);


figure
errorbar(1:16,mean(codeWords(c<-45,:)),std(codeWords(c<-45,:)))
hold on
errorbar(1.2:16.2,mean(codeWords(1:21,:)),std(codeWords(1:21,:)),'r')




figure
plot(1:16,codeWords(c<-44,:),'c')
hold on
plot(1:16,mean(codeWords(c<-44,:)),'b','linewidth',3)
plot(1:16,codeWords(1:25,:),'m')
plot(1:16,mean(codeWords(1:25,:)),'r','linewidth',3)

figure
colormap('gray')
c1=c;
c1(c<-44)=-100;
imagesc(reshape(c1,60,60))


figure
tr=44
plot(codeWords(c<-tr,6)+codeWords(c<-tr,13)-codeWords(c<-tr,4)-codeWords(c<-tr,3),'b')
hold on
plot(codeWords(1:25,6)+codeWords(1:25,13)-codeWords(1:25,4)-codeWords(1:25,3),'r')

figure
plot(c,codeWords(:,6)+codeWords(:,13)-codeWords(:,4)-codeWords(:,3),'.k')


figure
bar(codeWords(d,:))

figure
bar(mean(codeWords(1:100,:)))


figure
colormap('gray')
plot((std(codeWords)))


figure
colormap('gray')
imagesc(sta(:,:,4))



a=reshape(codeWords,40,40,16);
figure
colormap('gray')
imagesc(a(:,:,5))
figure
plot(codeWords)



for i=1:4
    subplot(2,2,i)
    colormap('gray')
    imagesc(sta(:,:,i))
end


[a,b]=find(std(sta,0,3)>140)
for i=1:length(a)
    tt(i,1:5)=reshape(sta(a(i),b(i),:),5,1);
end

templ=(mean(tt)/length(spikes)-0.5)'*50;
plot(templ)

figure
imagesc(sta(130:179,65:114,3))

cnt=1;
for i=125:184
    for j=60:119
        k(cnt)=sub2ind([200,200],i,j);
        cnt=cnt+1;
    end
end
tmp2(k(end),3);
sta(174,109,3)

tt=tmp2(k,:);
for i=1:4
    subplot(2,2,i)
    colormap('gray')
    imagesc(reshape(tt(:,i),40,40))
end