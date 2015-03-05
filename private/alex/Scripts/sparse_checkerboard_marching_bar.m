cd('S:\user\alexandra\scripts')

%% Assign colors - DO NOT RUN this part!
load('S:\user\alexandra\stimfiles\_hlc2_\goodGroupsArray.mat')
while max(abs(sum(reshape(matr,100,152))))>6
    matr(find(matr))=1;
    for i=1:10
        for j=1:10
            a=find(matr(i,j,:)==1);
            m=randperm(length(a));
            if mod(length(a),2)
                b=rand(1)<0.5;                    
            else
                b=0;
            end                   
            m(floor(length(a)/2)+1+b:end)=[];
            matr(i,j,a(m))=-1;
        end
    end    
    max(abs(sum(reshape(matr,100,152))))
end
a=sum(reshape(matr,100,152));
hist(a,max(a)-min(a))
title(sum(a))

sum(matr,3)
%%%%%%---------- DO NOT OVERWRITE THE EXISTING FILE!! ----------%%%%%%
% save('S:\user\alexandra\stimfiles\_hlc2_\CB_sparse.mat','matr')


%% Create stimuli pictures

load('S:\user\alexandra\stimfiles\_hlc2_\CB_sparse.mat','matr')

imageSize=670; % pxls
checkSize=67; % pxls
checkColor=10; % 'black' in RGB
checkColor2=50; % 'white' in RGB
bgcolor=30; % 'background' in RGB
pathName='S:\user\alexandra\stimfiles\_quickControl_\';
cnt=4; % from which number to start picture numbering

checkColor=checkColor*257;%convertion to 0-65535 scheme
checkColor2=checkColor2*257;%convertion to 0-65535 scheme
bgcolor=bgcolor*257;

coord=1:67:67*10; % coordinates of the bottom left checker pixels

for i=1:152
    s=strcat('00',num2str(cnt,'%d'));%the number of the picture for matlab stim directory
    s=s(end-2:end);
    fileName=strcat(s,'_cbSparse_pattern_',int2str(i),'.png');
    cnt=cnt+1;
    alpha=zeros(imageSize,imageSize)+65535;%px;670 covers about 2010µm - fully transparent (background is visible)
    alpha=uint16(alpha);
    
    bg=zeros(imageSize,imageSize)+bgcolor;
    [a,b]=find(matr(:,:,i));
    for j=1:length(a)
        bg(coord(a(j)):coord(a(j))+66,coord(b(j)):coord(b(j))+66)=matr(a(j),b(j),i);
    end
    bg(bg==-1)=checkColor;
    bg(bg==1)=checkColor2;

    im=cat(3,bg,bg,bg);
    im=uint16(im);
    imwrite(im, strcat(pathName,fileName),'Alpha',alpha);
end

%% Create stimuli protocols .stim.txt with 8 patterns per protocol

pathName='S:\user\alexandra\stimfiles\_quickControl_\';
cd(pathName);
cnt=3; % add to match the picture name number in the folder
for i=1:8:152
    fid=fopen([pathName,'cbSparse_pat_',int2str(i),'_',int2str(i+7),'.stim.txt'],'w');
    fprintf(fid,strcat('Parameters:\r\nID\tduration\tposx\tpoy\r\n-1\t500\t0\t0'));
    for j=i:i+7
        fprintf(fid,['\r\n',int2str(j+cnt),'\t2000\t0\t0\r\n']);
        fprintf(fid,'-1\t3000\t0\t0');
    end
    fclose(fid);
end




%% Create stimuli pictures: bars
pathName='S:\user\alexandra\stimfiles\_quickControl_\';

barShortSide=200; % pxls
imageSize=2000; % pxls
stimcolorBlack=10; % 'black' in RGB
stimcolorWhite=50; % 'white' in RGB
bgcolor=30;
pathName='S:\user\alexandra\stimfiles\_quickControl_\';
stimcolorBlack=stimcolorBlack*257;%convertion to 0-65535 scheme
stimcolorWhite=stimcolorWhite*257;%convertion to 0-65535 scheme
bgcolor=bgcolor*257;

alpha=zeros(imageSize,imageSize)+65535;%px;670 covers about 2010µm - fully transparent (background is visible)
alpha=uint16(alpha);
fileName='002_barBlackVert.png';
bg=zeros(imageSize,imageSize)+stimcolorBlack;
bg(:,(imageSize/2-barShortSide/2+1):(imageSize/2+barShortSide/2))=bgcolor;
im=cat(3,bg,bg,bg);
im=uint16(im);
imwrite(im, strcat(pathName,fileName),'Alpha',alpha);
fileName='003_barWhiteVert.png';
bg=zeros(imageSize,imageSize)+stimcolorWhite;
bg(:,(imageSize/2-barShortSide/2+1):(imageSize/2+barShortSide/2))=bgcolor;
im=cat(3,bg,bg,bg);
im=uint16(im);
imwrite(im, strcat(pathName,fileName),'Alpha',alpha);


%% Create stimuli protocols .stim.txt with 8 patterns per protocol: bars
pathName='S:\user\alexandra\stimfiles\_quickControl_\';
cd(pathName);

% order: black vertical (1:40, 1:8 - positions from left to right, then repetition),
% white vertical (41:90)
a=randperm(40);
a=[a; a+44];% total presentations of 4 bar types in random order
a=a(:);
coord=-1000:200:400;

for i=1:10:80
    fid=fopen([pathName,'bars_pat_',int2str(i),'_',int2str(i+9),'.stim.txt'],'w');
    fprintf(fid,strcat('Parameters:\r\nID\tduration\tposx\tpoy\r\n-1\t500\t0\t0'));
    for k=i:i+9
        j=a(k);
        shift=int2str(coord(mod(j-1,8)+1));
        if j<=40
            pic='2';
            tt(k)=coord(mod(j-1,8)+1);
            wt(k)=nan;
        else
            pic='3';
            tt(k)=nan;
            wt(k)=coord(mod(j-1,8)+1);
        end
        fprintf(fid,['\r\n',pic,'\t2000\t',shift,'\t0\r\n']);
        fprintf(fid,'-1\t8000\t0\t0');
    end
    fclose(fid);
end



%% create batch files

% basic set
cd('S:\user\alexandra\stimfiles\_quickControl_\');

fid_batch=fopen('basic8.batch.txt','w');
fprintf(fid_batch,'nd\tnone\t0\r\n');

for k=1:8
    fprintf(fid_batch,'HCStefano_1\tnone\t0\r\n');
    fprintf(fid_batch,'LCStefano_1\tnone\t0\r\n');
    fprintf(fid_batch,'quick\tnone\t0\r\n');
    fprintf(fid_batch,'chirp\tnone\t0\r\n');
    fprintf(fid_batch,'chirp\tnone\t0\r\n');
    fprintf(fid_batch,'chirp\tnone\t0\r\n');
end

fprintf(fid_batch,'##BeginStimParam##\r\nCounter\r\n');
for i=1:8
    fprintf(fid_batch,'\r\n');
    fprintf(fid_batch,'\r\n');
    fprintf(fid_batch,'\r\n');
    fprintf(fid_batch,'30\t30\t500\r\n');
    fprintf(fid_batch,'30\t30\t500\r\n');
    fprintf(fid_batch,'30\t30\t500\r\n');
end
fprintf(fid_batch,'##BeginVoltParam##');
for i=1:49
    fprintf(fid_batch,'\r\n');
end
fprintf(fid_batch,'##BeginCounterSettings##\r\n');
fclose(fid_batch);



% full set
cd('S:\user\alexandra\stimfiles\_quickControl_\');

fid_batch=fopen('full.batch.txt','w');
fprintf(fid_batch,'nd\tnone\t0\r\n');

%basic 2 times
for k=1:2
    fprintf(fid_batch,'HCStefano_1\tnone\t0\r\n');
    fprintf(fid_batch,'LCStefano_1\tnone\t0\r\n');
    fprintf(fid_batch,'quick\tnone\t0\r\n');
    fprintf(fid_batch,'chirp\tnone\t0\r\n');
    fprintf(fid_batch,'chirp\tnone\t0\r\n');
    fprintf(fid_batch,'chirp\tnone\t0\r\n');
end
%center 1 time
for k=1:8:152
    fprintf(fid_batch,['cbSparse_pat_',int2str(k),'_',int2str(k+7),'\tnone\t0\r\n'])
end
%basic 1 time
fprintf(fid_batch,'HCStefano_1\tnone\t0\r\n');
fprintf(fid_batch,'LCStefano_1\tnone\t0\r\n');
fprintf(fid_batch,'quick\tnone\t0\r\n');
fprintf(fid_batch,'chirp\tnone\t0\r\n');
fprintf(fid_batch,'chirp\tnone\t0\r\n');
fprintf(fid_batch,'chirp\tnone\t0\r\n');
%surround 1 time
for k=1:10:80
    fprintf(fid_batch,['bars_pat_',int2str(k),'_',int2str(k+9),'\tnone\t0\r\n'])
end
%basic 1 time
fprintf(fid_batch,'HCStefano_1\tnone\t0\r\n');
fprintf(fid_batch,'LCStefano_1\tnone\t0\r\n');
fprintf(fid_batch,'quick\tnone\t0\r\n');
fprintf(fid_batch,'chirp\tnone\t0\r\n');
fprintf(fid_batch,'chirp\tnone\t0\r\n');
fprintf(fid_batch,'chirp\tnone\t0\r\n');
%binary checkers 1 time
 fprintf(fid_batch,'check_bin\r\n');
%basic 6 times
for k=1:6
    fprintf(fid_batch,'HCStefano_1\tnone\t0\r\n');
    fprintf(fid_batch,'LCStefano_1\tnone\t0\r\n');
    fprintf(fid_batch,'quick\tnone\t0\r\n');
    fprintf(fid_batch,'chirp\tnone\t0\r\n');
    fprintf(fid_batch,'chirp\tnone\t0\r\n');
    fprintf(fid_batch,'chirp\tnone\t0\r\n');
end


fprintf(fid_batch,'##BeginStimParam##\r\nCounter\r\n');
for i=1:2
    fprintf(fid_batch,'\r\n');
    fprintf(fid_batch,'\r\n');
    fprintf(fid_batch,'\r\n');
    fprintf(fid_batch,'30\t30\t500\r\n');
    fprintf(fid_batch,'30\t30\t500\r\n');
    fprintf(fid_batch,'30\t30\t500\r\n');
end
for k=1:19
    fprintf(fid_batch,'\r\n');
end
fprintf(fid_batch,'\r\n');
fprintf(fid_batch,'\r\n');
fprintf(fid_batch,'\r\n');
fprintf(fid_batch,'30\t30\t500\r\n');
fprintf(fid_batch,'30\t30\t500\r\n');
fprintf(fid_batch,'30\t30\t500\r\n');
for k=1:8
    fprintf(fid_batch,'\r\n');
end
fprintf(fid_batch,'\r\n');
fprintf(fid_batch,'\r\n');
fprintf(fid_batch,'\r\n');
fprintf(fid_batch,'30\t30\t500\r\n');
fprintf(fid_batch,'30\t30\t500\r\n');
fprintf(fid_batch,'30\t30\t500\r\n');
fprintf(fid_batch,'\r\n');
for i=1:6
    fprintf(fid_batch,'\r\n');
    fprintf(fid_batch,'\r\n');
    fprintf(fid_batch,'\r\n');
    fprintf(fid_batch,'30\t30\t500\r\n');
    fprintf(fid_batch,'30\t30\t500\r\n');
    fprintf(fid_batch,'30\t30\t500\r\n');
end



fprintf(fid_batch,'##BeginVoltParam##');
for i=1:89
    fprintf(fid_batch,'\r\n');
end

fprintf(fid_batch,'##BeginCounterSettings##\r\n');
fclose(fid_batch);

