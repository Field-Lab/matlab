for i=1:10
    FileName=['C:\pawel\nauka\analiza\retina\2012-09-27-4\files\proby2\e12delay' num2str(i+1033) 'Neuron.tif'];
    RGB = imread(FileName);
    [X,map] = rgb2ind(RGB,128);
    OutputFileName='C:\pawel\nauka\analiza\retina\2012-09-27-4\files\proby2\e12proba.gif';
    if i == 1;
        imwrite(X,map,OutputFileName,'gif', 'Loopcount',inf);
    else
        imwrite(X,map,OutputFileName,'gif','WriteMode','append','DelayTime',0.1);
    end
end