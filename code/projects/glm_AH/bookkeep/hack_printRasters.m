function hack_printRasters(organizedspikes,SPars,pdf_name)


testseconds     = SPars.nsec_o;
plotseconds     = 10;
figs            = floor (testseconds / plotseconds) ;
RasterBlocks = SPars.StaticBlocks;
ytick = 0:10:length(RasterBlocks);      
clf; 
MS = 6;


spikes = 0;
for i_fit = SPars.NovelBlocks
    spikes = spikes + length(organizedspikes.block.t_sp{i_fit});
end

sprate_hz = spikes / ( length(SPars.NovelBlocks) * SPars.nsec_e);


clf

subplot((figs+1), 1, (1))
c = 0;
c = c+2;
text(-.1,1-0.1*c,sprintf('%s Cell %d',organizedspikes.cell_type, organizedspikes.cell_id))
c = c+2;
text(-.1,1-0.1*c,sprintf('Spike Rate of %1.3e Hz',sprate_hz))
c = c+2;
text(-.1,1-0.1*c,sprintf('%s',datestr(clock) ))
axis off

for i_fig = 1:figs
    subplot((figs+1), 1,(1+i_fig));        
    sec_start = 0 + (i_fig-1) *plotseconds;
    sec_end   = i_fig *plotseconds;
    xtick = sec_start:sec_end;

	for i_rep = 1:length(RasterBlocks);
        i_blk = RasterBlocks(i_rep);
        tsp0  = organizedspikes.block.t_sp_withinblock{i_blk}; 
        tsp0  = tsp0(tsp0 < sec_end);
        tsp   = tsp0(tsp0 > sec_start);
        yval  = i_rep * ones(size(tsp));
        plot(tsp,yval, 'r.','markersize',MS);
        if i_rep ==1
                hold on;
                set(gca, 'fontsize', 10);
                xlabel('Seconds');
                title('Recorded Spikes');
                ylabel('Repitiion Number');
                set(gca, 'xtick', xtick);
                set(gca, 'ytick', ytick);
                
                xlim([sec_start , sec_end])
                ylim([0, length(RasterBlocks)])
        end
	end
        
end
orient landscape

eval(sprintf('print -dpdf %s.pdf',pdf_name))
end


