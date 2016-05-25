[DataTraces,ArtifactDataTraces,Channels]=NS_ReadPreprocessedData(DataPath,ArtifactDataPath,0,45,118,1500,0);
el=11;
s=reshape(DataTraces(:,el,:),50,600);
ss=size(s);
figure(1)
%break
t=[1:600]/20;
FS=20
FigurePath='D:\Home\Pawel\analysis\slices\2013\2010-09-14-0\data002_Matlab\rozne\Anim_n228_stim64_m121';
for i=1:ss(1)
    plot(t,s(i,:));
    axis([0 30 -150 60]);
    grid on
    refresh
    h=text(26,50,['rep. ' num2str(i)])
    set(h,'FontSize',FS);
    xlabel('Time [ms]');
    h=gca
    set(h,'FontSize',FS);
    
    h=gcf;
    frame = getframe(h);
    im = frame2im(frame);
    [imind,map] = rgb2ind(im,256);
    if i == 1
        %imwrite(imind,map,FigurePath,'gif', 'Loopcount',1);
    else
        %imwrite(imind,map,FigurePath,'gif','WriteMode','append','DelayTime',0.2);
    end
    
    %imwrite(imind,map,[FigurePath '_f' num2str(i)],'gif');
    pause(0.2);
end