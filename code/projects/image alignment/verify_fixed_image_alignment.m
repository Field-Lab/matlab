% verify the aligment of cell bodies in fixed image and alive images
%
% plots image of living retina with cell body locations found in alive images (blue) and fixed images (red)
% the two versions of each cell are connected with a green line to verify they are the same cell


% PARAMETERS
fig_1 = 1;
fig_2 = 102;

alive_alignment = a1;
alive_image = a1;
alive_tform = TA1;



% plot points on fixed images and on alive alignment image separately
figure(fig_1);clf;
ss = 4;
subplot(121);imagesc(stacks(ss).max);axis equal;hold on;
plot(stacks(ss).ip(:,1),stacks(ss).ip(:,2),'.r');
subplot(122);imagesc(alive_alignment);axis equal;hold on;
plot(stacks(ss).bp(:,1),stacks(ss).bp(:,2),'.b')



% plot fixed and alive points on the same alive image

% transform alive image to array coordinates
[im_array,im_xdata,im_ydata] = imtransform(alive_image,alive_tform,'xdata',1100*[-1 1],'ydata',500*[-1 1]);
% plot it
figure(fig_2);clf;
imagesc(im_array,'xdata',im_xdata,'ydata',im_ydata);hold on;

for ss = 1:length(stacks)

    % get cell body locations from each fixed image
    cells_ip = tforminv(stacks(ss).tforminv,stacks(ss).ip);
    cells_bp = tformfwd(TA1,stacks(ss).bp);
    plot(cells_ip(:,1),cells_ip(:,2),'.r',cells_bp(:,1),cells_bp(:,2),'.b')

    % draw lines connecing them
    for cc=1:size(cells_ip,1)
        plot([cells_ip(cc,1) cells_bp(cc,1)],[cells_ip(cc,2) cells_bp(cc,2)],'g')
    end


    % compute rectangle of the image in array coordinates
    wd = getfield(imfinfo(stacks(ss).im{1}),'Width');
    ht = getfield(imfinfo(stacks(ss).im{1}),'Height');
    array_rect = tforminv(stacks(ss).tforminv,[1 1;1 ht; wd ht; wd 1; 1 1]);

    % plot it
    plot(array_rect(:,1),array_rect(:,2),'color',[.5 .5 0])

    % plot name of this
    text(mean(array_rect(1:4,1)),mean(array_rect(1:4,2)),sprintf('stacks(%d)',ss),'Color',[.5 .5 0],'FontSize',12,...
        'HorizontalAlignment','Center','VerticalAlignment','Bottom')

end

% add electrodes
if 1
    ep = datarun.ei.position;
    plot(ep(:,1),ep(:,2),'.','Color',[1 1 0])
    for ee=1:size(ep,1)
        %text(ep(ee,1),ep(ee,2),num2str(ee),'Color',[.5 .5 0],'FontSize',10,...
        %    'HorizontalAlignment','Center','VerticalAlignment','Bottom')
    end
end

% ensure everything is visible
axis image;axis equal
