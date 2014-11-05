function adjust_cones(frame, hTemplate)

global datarun cones myCells stim
persistent myPlot


tmp = get(hTemplate(1), 'SelectedObject');

stim=[];
stixel=datarun.stimulus.stixel_height;
stimarea=1;
stim.stimarea=stimarea;



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

if ishandle(myPlot)
    delete(myPlot)
end
myPlot=subplot('position',[0.6 0.2 0.3 0.3]);
set(gca,'DataAspectRatio',[1 1 1])
hold on
mycolors='rbkgmcykkkkkk';

pp=1;
col=mycolors(1);


if tmp==hTemplate(2) % gauss
    rf_size=[5 5];
    cone_templ=zeros(5,5,25);
    exactCorrection=zeros(25,2);
    cnt=1;
    for i=-2:2
        for j=-2:2
            bw_kern = make_gaussian('dim',2,'x_size',rf_size(2),'y_size',rf_size(1),'normalize','sum',...
                'center_radius',0.65,'center',[3+i/3,3+j/3], 'effective_radius',2);
            cone_templ(:,:,cnt) = full(bw_kern);
            exactCorrection(cnt,:)=[i/3,j/3];
            cnt=cnt+1;
        end
    end
    
    if min(all_cones(3,3,:))<0 % OFF cells
        all_cones=-all_cones; % make ON RFs
    end

    cnt=1; pp=1;
    for i=1:length(myCones)
        
        tmp=all_cones(:,:,i);
        %         figure
        %         imagesc(tmp)
        subCone=tmp(2:4,2:4);
        maxInd=max(subCone(:));
        tmp(tmp>maxInd)=0;
        maxCor=zeros(1,25);
        for j=1:25
            A = tmp .* cone_templ(:,:,j);
            maxCor(j) = sum(A(:));
        end
        
        [myweight, templateNumber]=max(maxCor);
        
        myCorrection=exactCorrection(templateNumber,:);
               
        
        a(1)=myCones(i,1)+myCorrection(1);
        a(2)=myCones(i,2)+myCorrection(2);
        
        stim.weight(i)=myweight;
        stim.templateNumber(i)=templateNumber;
        
        stim.coord(i,1)=a(1);
        stim.coord(i,2)=a(2);
        
        %     plot(a(2),a(1),'*r')
        plot([a(1)-stimarea a(1)+stimarea],[a(2)-stimarea a(2)-stimarea],col)
        
        plot([a(1)-stimarea a(1)+stimarea],[a(2)+stimarea a(2)+stimarea],col)
        
        plot([a(1)-stimarea a(1)-stimarea],[a(2)-stimarea a(2)+stimarea],col)
        
        plot([a(1)+stimarea a(1)+stimarea],[a(2)-stimarea a(2)+stimarea],col)
        
        if cnt>=conesPerCell(pp)
            pp=pp+1;
            col=mycolors(pp);
            cnt=1;
        else
            cnt=cnt+1;
        end
    end
 
    % check for shared cones
    w=squareform(pdist(stim.coord));
    w(w==0)=max(w(:));
    p=min(w);
    if ~isempty(find(p<1, 1)) % doubled cone
        p=find(p<1);
        i=1;
        while length(p)>1
            q=find(w(:,p(1))==min(w(:,p(1))));
            [~,n]=max(stim.weight([p(1) q]));
            if n==1
                y=[p(1),q];
            else y=[q,p(1)];
            end
            stim.coord(y(2),:)=stim.coord(y(1),:);
            stim.templateNumber(y(2))=stim.templateNumber(y(1));
            p(p==q)=[];
            i=i+1;
        end
    end
    
    axis([0 datarun.stimulus.field_width 0 datarun.stimulus.field_height])
    
else
    
    cone_templ=mean(all_cones,3);
    cone_templ=cone_templ/abs(cone_templ(3,3));
    
    
    cone_templ=imresize(cone_templ,stixel);
    
    
    all_resized=imresize(all_cones,stixel,'method','nearest');
    
    midStixel=size(all_resized,1)/2-0.5;
    
    
    
    
    for i=1:length(myCones)
        
        
        tmp=all_resized(:,:,i);
        tmp(tmp<all_cones(3,3,i))=0;
        
        A = conv2(tmp, cone_templ);
        [~, imax] = max(abs(A(:)));
        [ypeak, xpeak] = ind2sub(size(A),imax(1));
        a=[ypeak,xpeak]-midStixel;
        
        a(2)=myCones(i,1)*stixel - midStixel -2 + a(2);
        a(1)=myCones(i,2)*stixel - midStixel -2 + a(1);
        
        stim.coord(i,1)=a(2);
        stim.coord(i,2)=a(1);
        
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
    axis([0 stixel*datarun.stimulus.field_width 0 stixel*datarun.stimulus.field_height])
    
end
axis ij

stim.cone_template=cone_templ;

