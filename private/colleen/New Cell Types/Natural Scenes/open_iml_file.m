% open natural scenes image iml file
close all
for i = 23
    if i <10
        f1 = fopen(['imk0000', num2str(i),'.iml'], 'rb', 'ieee-be');
    else
                f1 = fopen(['imk000', num2str(i),'.iml'], 'rb', 'ieee-be');

    end
    
w = 1536; h = 1024;
buf = fread(f1, [w, h], 'uint16');
figure;
colormap(gray);

buf2 = buf(end-599:end-299, 1:300);
 imagesc(buf2');axis equal; axis off
end


