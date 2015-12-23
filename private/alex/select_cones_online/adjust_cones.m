function adjust_cones(frame)

global datarun cones myCells
persistent myPlot

stixel=datarun.stimulus.stixel_height;

myCones=cell2mat(cones');
all_cones=zeros(5,5,length(cell2mat(cones')));
cnt=1;
for i=1:length(myCells)
    datInd=find(datarun.cell_ids==myCells(i));
    sta=squeeze(datarun.stas.stas{datInd});
    sta=sta(:,:,frame); % frame
    
    for j=1:length(cones{i})
        all_cones(:,:,cnt)=sta(cones{i}(j,2)-2:cones{i}(j,2)+2,cones{i}(j,1)-2:cones{i}(j,1)+2);
        cnt=cnt+1;
    end
    conesPerCell(i)=size(cones{i},1);
end


cone_templ=mean(all_cones,3);
cone_templ=cone_templ/abs(cone_templ(3,3));
cone_templ=imresize(cone_templ,stixel);

all_resized=imresize(all_cones,stixel,'method','nearest');

midStixel=size(all_resized,1)/2-0.5;

% bb=imresize(sta,3,'nearest');
% figure
% colormap gray
% imagesc(bb)
% hold on

stimarea=4;

if ishandle(myPlot)
    delete(myPlot)
end
myPlot=subplot('position',[0.6 0.2 0.3 0.3]);
set(gca,'DataAspectRatio',[1 1 1])
hold on
mycolors='rbkgmcykkkkkk';

cnt=1;pp=1;
col=mycolors(1);
for i=1:length(myCones)
    
    
    
    tmp=all_resized(:,:,i);
    tmp(tmp<all_cones(3,3,i))=0;

    A = conv2(tmp, cone_templ);
    [~, imax] = max(abs(A(:)));
    [ypeak, xpeak] = ind2sub(size(A),imax(1));
    a=[ypeak,xpeak]-midStixel;
    
    a(2)=myCones(i,1)*stixel - midStixel -2 + a(2);    
    a(1)=myCones(i,2)*stixel - midStixel -2 + a(1);  

%     plot(a(2),a(1),'*r')
    plot([a(2)-stimarea a(2)+stimarea],[a(1)-stimarea a(1)-stimarea],col)

    plot([a(2)-stimarea a(2)+stimarea],[a(1)+stimarea a(1)+stimarea],col)
    
    plot([a(2)-stimarea a(2)-stimarea],[a(1)-stimarea a(1)+stimarea],col)
    
    plot([a(2)+stimarea a(2)+stimarea],[a(1)-stimarea a(1)+stimarea],col)
   
    if cnt>=conesPerCell(pp)
        pp=pp+1;
        col=mycolors(pp);
        cnt=1;
    else
        cnt=cnt+1;
    end
end
axis ij
axis([0 stixel*datarun.stimulus.field_width 0 stixel*datarun.stimulus.field_height])