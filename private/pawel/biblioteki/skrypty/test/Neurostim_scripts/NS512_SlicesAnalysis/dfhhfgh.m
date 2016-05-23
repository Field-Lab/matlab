for i=2:25
    RGB1 = imread(['D:\Home\Pawel\analysis\slices\2013\2013-12-15-0\SpikeAnalysis\data004\Amplitude' num2str(i) '.tif']);
    RGB2 = imread(['D:\Home\Pawel\analysis\slices\2013\2013-12-15-0\SpikeAnalysis\data008\Amplitude' num2str(i) '.tif']);
    RGB=RGB1;
    RGB(541:1080,:,:)=RGB2(1:540,:,:);
    IMWRITE(RGB,['D:\Home\Pawel\analysis\slices\2013\2013-12-15-0\SpikeAnalysis\data004008\Amplitude' num2str(i) '.tif'],'tiff');
end