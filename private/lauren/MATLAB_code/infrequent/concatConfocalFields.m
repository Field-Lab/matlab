
cd('/jacob/snle/data/2009-06-02-0/confocal/2009-06-12/')

imgString = cell(3,4);
imgString{1,4} = '2009-06-02-0_Series127';
imgString{2,4} = '2009-06-02-0_Series128';
imgString{3,4} = '2009-06-02-0_Series129';
imgString{3,3} = '2009-06-02-0_Series130';
imgString{2,3} = '2009-06-02-0_Series131';
imgString{1,3} = '2009-06-02-0_Series132';
imgString{1,2} = '2009-06-02-0_Series133';
imgString{2,2} = '2009-06-02-0_Series134';
imgString{3,2} = '2009-06-02-0_Series135';
imgString{3,1} = '2009-06-02-0_Series136';
imgString{2,1} = '2009-06-02-0_Series137';

% imgAllChan = cell(30,1);
% for i = 1:30
%     imgAllChan{i} = 255*ones(3072, 4096, 3, 'uint8');
% end


% img1 = imread('2009-06-02-0_Series127_z020_ch00.tif')/255;
% img2 = imread('2009-06-02-0_Series128_z020_ch00.tif')/255;
% img3 = imread('2009-06-02-0_Series129_z020_ch00.tif')/255;




for i = 1:30 %z stacks
    disp(['working on ' num2str(i)])
    imgAllChan = 255*ones(3072, 4096, 3, 'uint8');
    for j = 1:3 %rows
        for k = 1:4 %columns
            for m = 1:3 %channels
                if ~(j==1 && k == 1) %empty spot
                    if i < 11
                        tempImg = imread([imgString{j,k} '_z00' num2str(i-1) '_ch0' num2str(m-1) '.tif']);
                    else
                        tempImg = imread([imgString{j,k} '_z0' num2str(i-1) '_ch0' num2str(m-1) '.tif']);
                    end
                    imgAllChan(1024*(j-1)+1:1024*j, 1024*(k-1)+1:1024*k, m) = tempImg;
                end
            end
        end
    end
    if i < 11
        imwrite(imgAllChan, ['2009-06-02-0_Series127-137_Merged_z00' num2str(i-1)], 'tiff')
    else
        imwrite(imgAllChan, ['2009-06-02-0_Series127-137_Merged_z0' num2str(i-1)], 'tiff')
    end
end
