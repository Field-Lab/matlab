
% get data
if 1
    % good matches
    switch 1
        case 1 % 2007-09-18-4

            good_matches=[9 11 14 15 39 41 44 ];
            % not: 17 65
    end


    all_corr = cell(1,length(good_matches));
    cell_ids = cell(1,length(good_matches));
    
    for gg = 1:length(good_matches)
        [all_corr{gg},cell_ids{gg}] = find_matches(datarun,axons,good_matches(gg));
        if good_matches(gg) == 15; all_corr{gg} = all_corr{gg}([1 3:end]); cell_ids{gg} = cell_ids{gg}([1 3:end]); disp('here'); end
    end
end

figure(101);clf;hold on
for gg=1:length(good_matches)
    temp=sort(all_corr{gg},'descend');
    plot(temp(1),temp(2),'r.')
end
plot([0 1],[0 1],'k')
