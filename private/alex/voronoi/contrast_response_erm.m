function contrast_response_erm(inputs, spikes, nbins_cone1, nbins_cone2, center_cones, map, cones, comb, datarunID)

% contrast_bin_size in %: conditioning on the second cone
% inputs are filtered!

n_cones = size(inputs,1);
spikes(spikes>size(inputs,2)) = []; 
% get spike rate
spikes_tmp = spikes;
spike_rate=zeros(size(inputs,2),1);
while ~isempty(spikes_tmp)
    [~, ia, ~] = unique(spikes_tmp);
    spike_rate(spikes_tmp(ia))=spike_rate(spikes_tmp(ia))+1;
    spikes_tmp(ia)=[];
end
clear spikes_tmp

size_cone1_bin = floor(size(inputs,2)/nbins_cone1);
size_cone2_bin = floor(size(inputs,2)/nbins_cone2);

% on cell
% colors = [0.5 0.5 0.5; 1 0.5 0; 0.5 0 1;  1 1 0;0 1 0;  0 0 1;  1 0 0; 0 0.6 0; 1 0 1; 0 1 1; 0.6 0.65 0; ...
%     0.8 0.8 0; 0.8 0 0.8; 0 0.8 0.8; 0 0.3 0.8; 1 1 0.3; 0.9 0.3 0.8; 0.1 0.1 0.3; 0.8 0.8 0.8; 0.1 0.1 0.1];
% off cell
colors = [0.5 0.5 0.5;   1 1 0;0 1 0;  0 0 1;  1 0 0; 1 0 1; 0 1 1; 0.6 0.65 0];
    
for cone1 = 1:length(center_cones)
    
    figure(cone1)
    set(gcf, 'position', [154         369        1208         586])
    
    tmp = sort(inputs(cone1,:));    
    cone1_contr = tmp(1:size_cone1_bin:end);
    
    for cone2 = 1:size(inputs,1)
        if cone1 ~= cone2   
            cond_rate = zeros(nbins_cone2-1, nbins_cone1-1);
            if cone2<=length(center_cones)
                linewid = 3;
                cols = colors(cone2,:);
            else
                cols = [0 0 0];
                linewid=1;
            end
            % second cone: bin the inputs  
            tmp = sort(inputs(cone2,:));
            cone2_contr = tmp(1:size_cone2_bin:end);
             
             % go through contrast bins on cone2 and calculate the response
             % of cone 1
             cc=1;
             for i=[1,nbins_cone2]
                 % find instances when cone2 had certain contrast
                 cone2_valid = find(inputs(cone2,:)>cone2_contr(i) & inputs(cone2,:)<=cone2_contr(i+1));
                 
                 % go through contrast range on cone1
                 for j = 1:nbins_cone1
                     cone1_valid = find(inputs(cone1,:)>cone1_contr(j) & inputs(cone1,:)<=cone1_contr(j+1));
                     cone12_valid = intersect(cone1_valid, cone2_valid);
                     cond_rate(i,j) = mean(spike_rate(cone12_valid));
                 end  
                 subplot(2,3,cc)      
                 plot(cone1_contr(1:end-1)+diff(cone1_contr)/2,cond_rate(i,:)-mean(cond_rate(i,:)), 'linewidth', linewid, 'color', cols);
                 hold on
                 if cone1<length(center_cones)+1
                     title(['Ref cone ', int2str(center_cones(cone1))])
                 end
                 cc = cc+1;
             end  
            
        end
    end
    subplot(2,3,3)
    hold on
    leg = [];
    cnt = 1;
    for i=1:length(center_cones)
        if i~=cone1
%             if i<8
%                 marker='.';
%             else
%                 marker='x';
%             end
            t = plot(1:2, [i i], 'linewidth',3, 'color', colors(cnt,:));
            cnt = cnt+1;
            leg = [leg center_cones(i)];
        end
    end
    axis([0 3 0 i+1])
    legend(int2str(leg'))
    
    subplot(2,3,4)
    imagesc(comb)
    hold on
    for i=center_cones
        text(cones(i,2),cones(i,1), int2str(i), 'color', 'r')
    end
    coords = cones(center_cones,:);
    y = [min(coords(:,1))-10 max(coords(:,1))+10];
    x = [min(coords(:,2))-10 max(coords(:,2))+10];
    axis([x y])
%     
    
    % color cones
%     colored_cones = zeros([size(comb), 3]);  
    colored_cones = zeros([size(comb)]);  
    for cone2 = 1:length(center_cones)
        [a, b] = find(map==center_cones(cone2));        
        if cone1 ~= cone2            
            for j = 1:length(a)
                colored_cones(a(j),b(j), :) = colors(cone2,:);
            end
        else % ref cone
            for j = 1:length(a)
                colored_cones(a(j),b(j),:) = 1;
            end
        end
    end
    subplot(2,3,5)
    imagesc(colored_cones);
%     hold on
%     for cone2 = 1:length(center_cones)
%         if cone1 ~= cone2
%             text(cones(center_cones(cone2),2),cones(center_cones(cone2),1), int2str(center_cones(cone2)), 'color', [1 1 1])
%             text(cones(center_cones(cone2),2),cones(center_cones(cone2),1)+2, int2str(cone2), 'color', [1 1 1])
%         else
%             text(cones(center_cones(cone2),2),cones(center_cones(cone2),1), int2str(center_cones(cone2)), 'color', 'r')
%             text(cones(center_cones(cone2),2),cones(center_cones(cone2),1)+2, int2str(cone2), 'color', 'r')
%         end
%     end
    axis([x y])
    set(gca, 'dataaspectratio', [1 1 1])
    
    subplot(2,3,6)
    imagesc(comb)
    hold on

    
    drawnow
    path2save = ['/Users/alexth/Desktop/for_talk/', int2str(datarunID)];
    if ~isdir(path2save)
        mkdir(path2save);
    end
    saveas(gcf,[path2save,'/off_', int2str(center_cones(cone1)), '.tiff'])
end


