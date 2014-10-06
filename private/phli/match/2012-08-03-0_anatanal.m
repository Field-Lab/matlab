%% Just copied over from 2011-07-05-4, 2012-04-13-1; need to convert everything over.

%%
piece = '2012-08-03-0';
rig = 'A';
loadopts = struct('load_params', true, 'load_neurons', true, 'load_ei', true);
keep_vars = {'piece'; 'loadopts'; 'keep_vars'};


%%
d00 = load_data([piece '/streamed/data000/data000'], loadopts);
d00 = restore_stacks(d00);
% 
datarun = d00;
 
keep_vars = [keep_vars; 'd00'; 'datarun'];
 
d03.piece.array_id = 1504;
d03.piece.rig = 'A';


%% Setting up stacks; only run once then use restore_stacks
%datarun = init_stacks(datarun);
%datarun.stacks{2,1}  = build_stack([piece '/spot/low/stitch.tif']);
%datarun.stacks{2,2}  = build_stack([piece '/spot/hi/Stitch.tif']);
% 
% datarun = array_align(datarun, {2,1}, 'rig', rig);
% datarun = array_align(datarun, {2,2}, 'rig', rig);
% 
% datarun.stacks{10,1} = build_stack([piece '/confocal/2012-08-03-0_Brn3_rat/Brn3_rat_Series008_Merged_stack.tif:[1-7]']);
% datarun.stacks{10,2} = build_stack([piece '/confocal/2012-08-03-0_Brn3_rat/Brn3_rat_Series050_Merged_stack.tif:[1-10]']);
% datarun.stacks{10,3} = build_stack([piece '/confocal/2012-08-03-0_Brn3_rat/Brn3_rat_Series063_074_Merged_stack.tif:[1-16]']);
% 
% 
% save_stacks(datarun)

figure; imshow(get_slice(datarun.stacks{10,2}))


%% Reslice stack
reslicestack = {10,3};

% Pick points
datarun.stacks{reslicestack{:}} = load_slices(datarun.stacks{reslicestack{:}});
hreslicer = stack_reslicer(datarun.stacks{reslicestack{:}}, 'points_color', 'y', 'bgpoint_color', 'c', 'method', 'linear');
waitdata = waitfordata(hreslicer, 'points');
datarun.stacks{reslicestack{:}}.reslice_points = waitdata{1};
clear waitdata;


%%
save_stacks(datarun);


%% Save resliced image

resliced = getappdata(gcf, 'resliced');
[p,f] = fileparts(fullfile(datarun.stacks{reslicestack{:}}.basepath, datarun.stacks{reslicestack{:}}.paths{1}));
imwrite(resliced, sprintf('%s/%s%s.tif', p, f, '_resliced'));

% Now loop back up and add to stacks


%%
% stack_align
%%

% openall load_slices
% d03.stacks{11,6} = load_slices(d03.stacks{11,6})
% gui_points=[]
% stack_reslicer(d03.stacks{11,6}, 'points', gui_points, 'points_color', 'y');


%% Prep for plot physiology over image
datarun = get_lazyhex_ei_neighbors(datarun);
datarun = calc_csds(datarun, {1});


%% Plot ei over image

% cellid = d03.cell_types{1}.cell_ids(2);
% 
% plot_ei_scroll(d03, cellid, 'coordinates', 'stack', 'stack', {11,9}, 'max_scale', 0.75, 'scale', 1.5, 'pos_color', 'g');
% eiax = gca;
% set(eiax, 'Color', 'none');
% 
% stackax = axes;
% imshow(resliced(:,:,1));
% stackax = gca;
% 
% linkaxes([eiax stackax]);
% axes(eiax);
% 
% set(gcf, 'Units', 'normalized', 'Position', [0.5 0 0.5 1]);
% autozoom_ei(d03, cellid, 'stack', {11,3});
% axis fill;
% 


%% Plot current source densities
% celltype = 1;
% 
% ei_modification = @(ei,datarun,cellid)(ei2csd(ei, d03.ei.neighbor_struct));
% for i = 1:length(d03.cell_types{celltype}.cell_ids)
%     cellid = d03.cell_types{celltype}.cell_ids(i);
%     plot_ei_scroll(d03, cellid, 'modify_ei', ei_modification, 'max_scale', 0.75, 'pos_color', 'g', 'figure', i);
%     set(gca, 'Color', 'none');
% end

% Nice ON Parasols: 2 5 34
% Nice OFF Parasols: 5 9 26
% Nice ON Midgets: 6 11 16 36 51 66 71


%% Try to model fit EI or CSD
% First attempt: just two round 2D gaussians with gaussian timecourse


%% Compare EI and CSD
% celltype = 1;
% xlims = [10 25];
% 
% for i = 1:length(d03.cell_types{celltype}.cell_ids)
%     cellid = d03.cell_types{celltype}.cell_ids(i);
%     ei = get_ei(d03, cellid);
%     csd = ei2csd(ei, d03.ei.neighbor_struct);
% 
%     figure
%     sanesubplot(1,2,1);
%     plot(ei');
%     set(gca, 'XLim', xlims);
%     grid on;
%     
%     sanesubplot(1,2,2);
%     plot(csd');
%     set(gca, 'XLim', xlims);
%     grid on;
% end



%% Plot current source density over image

% cellid = d03.cell_types{1}.cell_ids(2);
% 
% ei_modification = @(ei,datarun,cellid)(ei2csd(ei, d03.ei.neighbor_struct));
% plot_ei_scroll(d03, cellid, 'coordinates', 'stack', 'stack', {11,9}, 'modify_ei', ei_modification, 'max_scale', 0.75, 'scale', 3, 'pos_color', 'g');
% eiax = gca;
% set(eiax, 'Color', 'none');
% 
% stackax = axes;
% imshow(resliced(:,:,1));
% stackax = gca;
% 
% linkaxes([eiax stackax]);
% axes(eiax);
% 
% set(gcf, 'Units', 'normalized', 'Position', [0.56 0 0.45 1]);
% autozoom_ei(d03, cellid, 'stack', {11,3}, 'padding_factor', 5, 'keep_aspect_ratio', false);
% %axis fill;
