FiguresPath='C:\pawel\nauka\analiza\slices\2010-09-14-0\analysis_2014_08_19\SpikesPropagation\';
AnimationPath='C:\pawel\nauka\analiza\slices\2010-09-14-0\analysis_2014_08_19\Anim1.gif';

for Frame=1:300
    FigureName=['13f' num2str(Frame) '.tif']
    a=imread([FiguresPath FigureName],'tif');
    [imind,map] = rgb2ind(a,256);    
    if Frame == 1
        imwrite(imind,map,AnimationPath,'gif', 'Loopcount',1);
    else
        imwrite(imind,map,AnimationPath,'gif','WriteMode','append','DelayTime',0);
    end
end
break
[imind,map] = rgb2ind(im,256);
if Frame == 1
    imwrite(imind,map,FigurePath,'gif', 'Loopcount',1);
else
    imwrite(imind,map,FigurePath,'gif','WriteMode','append','DelayTime',0);
end