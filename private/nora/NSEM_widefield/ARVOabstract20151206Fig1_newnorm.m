%% UAED TO MAKE FIG 2
clear
Corr_NS = [];
Loc = [];
Color_idx = [];
modu=[];
no_modu= [];
default_colors = get(gca,'ColorOrder');


%%
for cell_type = {'On Parasol', 'Off Parasol'}
    
    %%
    clear datarun NSEM_Corr locations stim_end
    stim_end = 20 - [240]/8+30;
    datarun{1} = load_data('/Volumes/Analysis/2015-09-23-0/data012-data013-data017-norefit/data017-from-data012_data013_data017/data017-from-data012_data013_data017');
    datarun{2} = load_data('/Volumes/Analysis/2015-09-23-0/data012-data013-data017-norefit/data012-from-data012_data013_data017/data012-from-data012_data013_data017');
    datarun{3} = load_data('/Volumes/Analysis/2015-09-23-0/data012-data013-data017-norefit/data013-from-data012_data013_data017/data013-from-data012_data013_data017');
    
    [NSEM_Corr,locations, prof1] = wide_field(datarun, cell_type, stim_end, 0);
    Corr_NS = [Corr_NS; NSEM_Corr(:,3)];  
    Loc = [Loc ; locations];
    Color_idx = [Color_idx; 1*ones(length(locations),1)];


    
    %%
    clear datarun NSEM_Corr locations stim_end
    stim_end = 20 -[280]/8+30;
    datarun{1} = load_data('/Volumes/Analysis/2015-10-06-0/data000-data015-norefit/data000-from-data000_data001_data002_data003_data004_data005_data006_data007_data008_data009_data010_data011_data012_data013_data014_data015/data000-from-data000_data001_data002_data003_data004_data005_data006_data007_data008_data009_data010_data011_data012_data013_data014_data015');
    datarun{2} = load_data('/Volumes/Analysis/2015-10-06-0/data000-data015-norefit/data004-from-data000_data001_data002_data003_data004_data005_data006_data007_data008_data009_data010_data011_data012_data013_data014_data015/data004-from-data000_data001_data002_data003_data004_data005_data006_data007_data008_data009_data010_data011_data012_data013_data014_data015');
    datarun{3} = load_data('/Volumes/Analysis/2015-10-06-0/data000-data015-norefit/data005-from-data000_data001_data002_data003_data004_data005_data006_data007_data008_data009_data010_data011_data012_data013_data014_data015/data005-from-data000_data001_data002_data003_data004_data005_data006_data007_data008_data009_data010_data011_data012_data013_data014_data015');
    
    [NSEM_Corr,locations, prof2] = wide_field_newnorm(datarun, cell_type, stim_end, 0);
    Corr_NS = [Corr_NS; NSEM_Corr(:,3)];
    Loc = [Loc ; locations];
    Color_idx = [Color_idx; 2*ones(length(locations),1)];


    
    %%
    
    clear datarun NSEM_Corr locations stim_end
    stim_end = 20 -[240, 320]/8+30;
    % datarun 1 = class, 2 = full_rep, 3 = left third, 4 = left two thirds
    datarun{1} = load_data('/Volumes/Analysis/2015-10-29-2/data045-data048/data048/data048');
    datarun{2} = load_data('/Volumes/Analysis/2015-10-29-2/data045-data048/data045/data045');
    datarun{3} = load_data('/Volumes/Analysis/2015-10-29-2/data045-data048/data046/data046');
    datarun{4} = load_data('/Volumes/Analysis/2015-10-29-2/data045-data048/data047/data047');
    [NSEM_Corr,locations, prof3] = wide_field_newnorm(datarun, cell_type, stim_end, 0);
    Corr_NS = [Corr_NS; NSEM_Corr(:,3); NSEM_Corr(:,4)];
    Loc = [Loc ; locations];
    Color_idx = [Color_idx; 3*ones(length(locations),1)];


    
%     %%
    clear datarun NSEM_Corr locations stim_end
    stim_end = 20 -[240, 320]/8+30;
    % datarun 1 = class, 2 = full_rep, 3 = left third, 4 = left two thirds
    datarun{1} = load_data('/Volumes/Analysis/2015-10-29-7/data008-data012/data012/data012');
    datarun{2} = load_data('/Volumes/Analysis/2015-10-29-7/data008-data012/data008/data008');
    datarun{3} = load_data('/Volumes/Analysis/2015-10-29-7/data008-data012/data009/data009');
    datarun{4} = load_data('/Volumes/Analysis/2015-10-29-7/data008-data012/data010/data010');
    [NSEM_Corr,locations, prof4] = wide_field_newnorm(datarun, cell_type, stim_end, 0);
    Corr_NS = [Corr_NS; NSEM_Corr(:,3); NSEM_Corr(:,4)];
    Loc = [Loc ; locations];
    Color_idx = [Color_idx; 4*ones(length(locations),1)];

    
    
end

%%
cutoff = 7;
figure;  hold on;
%fill([prof1(:,2)-(95/8); prof1(:,2)-(95/8)]+4,0.01+[prof1(:,1)/10000; zeros(200,1)], 0.2+[0.7 0.7 0.7],'EdgeColor', 0.2+[0.7 0.7 0.7])
x = prof1(:,2)-(95/8)+4;
offset = 0.028;
y = offset+prof1(:,1)/10000;
plot(x(x>-1.5),y(x>-1.5), 'Color', [0.7 0.7 0.7]-0.2)
plot(x(x>-1.5), zeros(length(x(x>-1.5)))+offset, 'Color', [0.7 0.7 0.7]-0.3);
plot(Loc(Loc>-cutoff), Corr_NS(Loc>-cutoff),'.', 'Color', 0.3*[1 1 1])
avg_rf = mean([prof1(:,1), prof2(:,1), -prof3(:,1),-prof4(:,1)],1);
hold on;
%plot(prof1(:,2)-(95/8),prof1(:,1)/4000, 'k', 'LineWidth', 1)
%plot(prof1(:,2)-(95/8),ones(200,1), 'k', 'LineWidth', 1)

plot([0 0], [-0.1 0.1], 'Color', 0.5*[1 1 1], 'LineWidth', 2)
% plot([-11 15], [0 0], 'k', 'LineWidth', 1)

[x, y, error, bincount, binedge] = curve_from_binning(Loc(Loc>-cutoff), Corr_NS(Loc>-cutoff), 'num_bins',8);
plot(x,y,'k', 'LineWidth', 2)


%ylim([-0.04 0.08])
set(gca, 'YTick', [0 0.03])
set(gca, 'XTick', [-5 0 5])
xlim([-7.5 7.5])
ylim([-0.001 0.039])

%% other not as good data
% try
% clear datarun NSEM_Corr locations stim_end
% stim_end = 20 -[260]/8+30;
% % datarun 1 = class, 2 = full_rep, 3 = left third, 4 = left two thirds
% datarun{1} = load_data('/Volumes/Analysis/2015-11-09-1/data009-data013/data009/data009');
% datarun{2} = load_data('/Volumes/Analysis/2015-11-09-1/data009-data013/data012/data012');
% datarun{3} = load_data('/Volumes/Analysis/2015-11-09-1/data009-data013/data013/data013');
% [NSEM_Corr,locations] = wide_field(datarun, cell_type, 0);
% datarun{1} = load_params(datarun{1});
% datarun{1}.default_sta_fits = 'vision';
% [r,m,s] = average_radius(datarun{1}, cell_type);
% Corr_NS = [Corr_NS; NSEM_Corr(:,3)];
% Loc = [Loc ; (locations(:,2)-stim_end)/(2*r)];
% Color_idx = [Color_idx; 5*ones(length(locations), 1)];
% figure(10);hold on; plot((locations(:,2)-stim_end)/(2*r), NSEM_Corr(:,3), '.', 'Color', default_colors(5,:))
% catch
%      disp('5');
% end

% clear datarun NSEM_Corr locations stim_end
% datarun{1} = load_data('/Volumes/Analysis/2015-08-17-2/data018-data022/data022-from-data018_data019_data020_data021_data022/data022-from-data018_data019_data020_data021_data022');
% datarun{2} = load_data('/Volumes/Analysis/2015-08-17-2/data018-data022/data020-from-data018_data019_data020_data021_data022/data020-from-data018_data019_data020_data021_data022');
% datarun{3} = load_data('/Volumes/Analysis/2015-08-17-2/data018-data022/data021-from-data018_data019_data020_data021_data022/data021-from-data018_data019_data020_data021_data022');
% datarun{4} = load_data('/Volumes/Analysis/2015-08-17-2/data018-data022/data018-from-data018_data019_data020_data021_data022/data018-from-data018_data019_data020_data021_data022');
% datarun{5} = load_data('/Volumes/Analysis/2015-08-17-2/data018-data022/data019-from-data018_data019_data020_data021_data022/data019-from-data018_data019_data020_data021_data022');
% [NSEM_Corr,locations] = wide_field(datarun, cell_type, 0);
%
%
% %% 30 um array, not many cells at all
% clear datarun NSEM_Corr locations stim_end
% datarun{1} = load_data('/Volumes/Analysis/2015-08-17-6/data022-24-repeats/data013-from-data013_data022_data023_data024/data013-from-data013_data022_data023_data024');
% datarun{2} = load_data('/Volumes/Analysis/2015-08-17-6/data022-24-repeats/data022-from-data013_data022_data023_data024/data022-from-data013_data022_data023_data024');
% datarun{3} = load_data('/Volumes/Analysis/2015-08-17-6/data022-24-repeats/data023-from-data013_data022_data023_data024/data023-from-data013_data022_data023_data024');
% datarun{4} = load_data('/Volumes/Analysis/2015-08-17-6/data022-24-repeats/data024-from-data013_data022_data023_data024/data024-from-data013_data022_data023_data024');
% [NSEM_Corr,locations] = wide_field(datarun, cell_type, 0);


