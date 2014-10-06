%[im,map] = rgb2ind(f.cdata,256,'nodither');
height=768;
width=1024;
frames=100;



im=zeros(height,width,1,frames);
frame0=ones(height,width,1,1)*200;
for i=1:frames
    frame=frame0;
    frame(i*2:i*2+30,:,1,1)=0;
    im(:,:,1,i)=frame;
    %imwrite(frame,['ble' num2str(i) '.gif'],'DelayTime',0,'LoopCount',inf);
    
    indX=round(rand*32)+1
    indY=round(rand*32)+1
    
    frame2=frame0;
    frame2(indX*30:(indX+1)*30,indY*30:(indY+1)*30)=0;
    imwrite(frame2,['mpg' num2str(i) '.jpg']);
end
%imwrite(im,'ble.gif','DelayTime',0,'LoopCount',inf);