function show_sta_IR(datarun,v,c,myCell)

persistent center_mark pixel_mark

a=floor(ginput(1));

for i=1:length(c)
    x=v(c{i}, 1); y=v(c{i}, 2);
    if sum(isinf(x))
        continue
    else
        IN = inpolygon(a(1),a(2),v(c{i}, 1),v(c{i}, 2));
        if IN
            break
        end
    end
end

a(1)=datarun.cones.centers(i,1);
a(2)=datarun.cones.centers(i,2);

if ishandle(center_mark)
    delete(center_mark)
end

center_mark=plot(a(1),a(2),'rx');

a=round(a);

if ishandle(pixel_mark)
    delete(pixel_mark)
end
pixel_mark=plot(a(1),a(2),'c+');

figure
hold on
% cone center STA
b=datarun.stas.stas{datarun.cell_ids==myCell};
if size(b,3)==1
    f=squeeze(b(a(2),a(1),1,:));
    
    plot(f)    
    f=squeeze(b(a(1)-1:a(1)+1,a(2)-1:a(2)+1,1,:));
    f=mean(reshape(f,9,size(b,4)));
    plot(f,'r')
    legend('cone center','3x3 cone')
else
    f=[squeeze(b(a(1),a(2),1,:))';squeeze(b(a(1),a(2),2,:))';squeeze(b(a(1),a(2),3,:))']';   
    colorlist='rgb';
    for i=1:3
        plot(f(:,i),colorlist(i))
    end
    
    f=squeeze(b(a(1)-1:a(1)+1,a(2)-1:a(2)+1,1,:));
    f=mean(reshape(f,9,size(b,4)));
    t=f';
    f=squeeze(b(a(1)-1:a(1)+1,a(2)-1:a(2)+1,2,:));
    f=mean(reshape(f,9,size(b,4)));
    t=[t f'];    
    f=squeeze(b(a(1)-1:a(1)+1,a(2)-1:a(2)+1,3,:));
    f=mean(reshape(f,9,size(b,4)));
    t=[t f'];
    
    for i=1:3
        plot(t(:,i),colorlist(i),'linewidth',2)
    end
    
    legend('center RED','center GREEN','center BLUE','3x3 RED','3x3 GREEN','3x3 BLUE')
    
end


line([1,size(b,4)],[0 0],'color','k')
