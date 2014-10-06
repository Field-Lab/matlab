clear all

datarun = load_data('2011-01-11-0/data030-lh/data030-lh');
%datarun = load_index(datarun);

datarun = load_params(datarun,'verbose',1);
%datarun = load_sta(datarun, 'save_rf', 1);
%datarun = load_ei(datarun, []);
datarun = load_neurons(datarun);

%cell_IDs{1} = [2 136 137 152 257 407 436 453 634 661 768 813 826 871]; %ON parasol
%cell_IDs{2} = [91 286 302 391 424 556 586 631 678 708 767 856 902]; %ON midget
%cell_IDs{3} = [18 106 138 288 392 451 531 573 617 751 753 888 905]; %OFF parasol
%cell_IDs{4} = [48 62 111 153 169 187 259 305 338 381 395 444 457 470 578 624 647 680 754 755 799 814 877 879 951]; %OFF midget
%cell_IDs{1} = datarun.cell_types{1}.cell_ids;

%all_IDs = [cell_IDs{1} cell_IDs{2} cell_IDs{3} cell_IDs{4}];
%datarun = get_sta_summaries(datarun, all_IDs);

%plot_rf_summaries(datarun, cell_IDs{1}, 'plot_fits', true)



%1 (ON P) 2 (OFF P) 3 (ON M) 4 (OFF M)

%%

for ii = 1:4
    acfs{ii} = [];

    %figure
    %hold on

    for jj = 1:length(datarun.cell_types{ii}.cell_ids)
        
        [time acf] = get_correlation(datarun, datarun.cell_types{ii}.cell_ids(jj));
        acf = acf/sqrt(sum(acf.^2)); %normalizes acf by rms

        %plot(1000*time, acf)

        
        acfs{ii} = [acfs{ii}; acf(end:-1:1)]; %reverse x-axis of plot
    end
    acfMeans{ii} = mean(acfs{ii},1);
    acfSDs{ii} = std(acfs{ii});
    %hold off
end

figure
hold on

colors = [1 0 0
          0 1 0
          1 0.5 0
          0 0 1];

for ii = 1:4     
    sdCloudY = [acfMeans{ii}+acfSDs{ii} acfMeans{ii}(end:-1:1)-acfSDs{ii}(end:-1:1)];
    sdCloudX = [time*1000 time(end:-1:1)*1000]; %time is in seconds, so convert to ms
    
    patch(sdCloudX, sdCloudY, colors(ii,:), 'edgeColor', 'none')
end

for ii = 1:4
    plot(time*1000, acfMeans{ii}, 'color', colors(ii,:), 'linewidth', 2)
end

xlabel('time difference (ms)')
ylabel('normalized correlation')




%% clusters

r=[];rr=[];rrr=[]; c=[];
for i = 1:4;
    %index = datarun.cell_types{i}.cell_ids;
    c = [c; repmat(colors(i,:), size(acfs{i},1), 1)];
    %     for ii=1:length(acfs{ii})
    r=[r acfs{i}(:,4:8)']; %one window of acf
    rr=[rr acfs{i}(:,8:end)']; %another window of acf
    rrr=[rrr acfs{i}(:,3:end)'];
    %
    %     end
end

[coeff,scores1] = princomp(r');
[coeff,scores2] = princomp(rr');
[coeff,scores3] = princomp(rrr');


%subplot(3,2,(j-1)*2+2)
figure
hold on

for ii=1:length(scores1(:,1))
    
    %plot(scores1(ii,1),scores1(ii,2),'.','color',c(ii,:),'MarkerSize',12)%12
    %plot(scores2(ii,1),scores2(ii,2),'.','color',c(ii,:),'MarkerSize',12)%12
    plot(scores3(ii,1),scores3(ii,2),'.','color',c(ii,:),'MarkerSize',12)%12

    %plot(scores5(ii,1),scores4(ii,1),'.','color',c(ii,:),'MarkerSize',5)%12
    %plot(scores2(ii,1),scores3(ii,1),'.','color',c(ii,:),'MarkerSize',5)%12
    %plot3(scores2(ii,1),scores2(ii,2),scores3(ii,1),'.','color',c(ii,:),'MarkerSize',5)%12
    hold on

end
set(gca,'PlotBoxAspectRatio',[1 1 1],'Xtick',[],'Ytick',[]);
%set(gca,'PlotBoxAspectRatio',[10 16 1],'Xtick',[],'Ytick',[]);
%[xrange yrange]=autoscale(scores3(:,1),scores2(:,1),'border', .15);
%set(gca,'XLim',xrange,'YLim',yrange);
box on

hold off






