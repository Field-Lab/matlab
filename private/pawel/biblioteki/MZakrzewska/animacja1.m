OutputFileName1='C:\pawel\nauka\dydaktyka\M_Zakrzewska\anim.gif';
figure(10);
for f=1:5
    %figure(f)
    %fgr=1.1^f
    plot(rand(1,100));
    h=gcf;
        
    frame = getframe(h);
    im = frame2im(frame);
    [imind,map] = rgb2ind(im,256);
    if f == 1
        imwrite(imind,map,OutputFileName1,'gif', 'Loopcount',1);
    else
        imwrite(imind,map,OutputFileName1,'gif','WriteMode','append','DelayTime',0.2);
    end
end