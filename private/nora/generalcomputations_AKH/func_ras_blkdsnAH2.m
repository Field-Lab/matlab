% AKHeitmans version of this

%%% Original Commenting
% plot rasters for a block-design experiment.
% n_blk blocks, odd blocks for the repeats; even blocks for non-repeats (aka long).
% 30 sec/blocks.
% based on ./chk_raster_blockdsn.m, extended for a) arbitrary block number; b) BW or NSEM stimulus.
% edoi, 2012-02-08.
%%%%%%%%%%%%%%%%%%%%%%%%%%%


function func_ras_blkdsnAH2(testrun,pdf_name)
%%

n_blk  = length(testrun.block.t_sp);
n_rep  = n_blk/2;
ax_rep = [0.5,n_rep+0.5];

clf;
c1 = 0;
c2 = 0;
sc = cell(2,1);
sc{1} = nan(n_rep,1);
sc{2} = nan(n_rep,1);
%sc_len = 30; % time length to count spikes
MS = 3;

ntb_o = testrun.no_trig_bl(1);
ntb_e = testrun.no_trig_bl(2);
%ntb_oe = ntb_o+ntb_e;

t_sp = testrun.t_sp;
t_trig = testrun.t_trig;

t_bin_wid = 50/1000;

t_bin1 = 0:t_bin_wid:ceil(testrun.no_trig_bl(1)*100/120);
t_bin2 = 0:t_bin_wid:ceil(testrun.no_trig_bl(2)*100/120);

t_itv1 = t_bin1(end)-1;
t_itv2 = t_bin2(end)-1;

PSTH = cell(2,1);
PSTH{1} = zeros(size(t_bin1));
PSTH{2} = zeros(size(t_bin2));
tot_n_sp = cell(2,1);


for k = 1:n_blk
   
   if rem(k,2) == 1
      c1 = c1+1;
      
      itv = [t_trig( (c1-1)*ntb_o + c2*ntb_e + 1 ), t_trig( c1*ntb_o + c2*ntb_e )];
      itv_idx = ( (t_sp > itv(1)) + (t_sp < itv(2)) == 2);
      itv_t_sp = t_sp(itv_idx);
      
      subplot(4,1.1,1)
      plot(itv_t_sp-t_trig( (c1-1)*ntb_o + c2*ntb_e + 1 ), c1*ones(size(itv_t_sp)),'k.','markersize',MS)
      sc{1}(c1) = sum(itv_t_sp-t_trig( (c1-1)*ntb_o + c2*ntb_e + 1 ) < t_itv1);
      tmp = histc(itv_t_sp-t_trig( (c1-1)*ntb_o + c2*ntb_e + 1 ),t_bin1);
      tmp = tmp(:)';
      PSTH{1} = PSTH{1}+ tmp;
      if c1 == 1
         hold on
         title(sprintf('%s (id%4d): Repeat blocks',testrun.cell_type,testrun.cell_ids))
         ylabel('Trial')
         set(gca,'ytick',[1,n_rep])
         ylim(ax_rep)
         xlim([0,t_itv1])
         axis ij
      end
   else
      c2 = c2+1;

      itv = [t_trig( c1*ntb_o + (c2-1)*ntb_e + 1 ), t_trig( c1*ntb_o + c2*ntb_e )];
      itv_idx = ( (t_sp > itv(1)) + (t_sp < itv(2)) == 2 );
      itv_t_sp = t_sp(itv_idx);

      
      subplot(4,1.1,1+2*1.1)
      plot(itv_t_sp-t_trig( c1*ntb_o + (c2-1)*ntb_e + 1 ), c2*ones(size(itv_t_sp)),'k.','markersize',MS)
      sc{2}(c2) = sum(itv_t_sp-t_trig(c1*ntb_o + (c2-1)*ntb_e + 1) < t_itv2);
      tmp = histc(itv_t_sp-t_trig(c1*ntb_o + (c2-1)*ntb_e + 1),t_bin2);
      PSTH{2} = PSTH{2}+tmp(:)';
      if c2 == 1
         hold on
         title('Random blocks')
         ylabel('Trial')
         set(gca,'ytick',[1,n_rep])
         ylim(ax_rep)
         xlim([0,t_itv2])
         axis ij
      end
   end
end

%%
%--
mx_y_psth = 0;
for i = 1:2
   tot_n_sp{i} = sum(PSTH{i});
   PSTH{i} = PSTH{i}/n_rep/t_bin_wid; % converting the unit to sps
   mx_y_psth = max([mx_y_psth,PSTH{i}]);
end

YL = [0,10*ceil(mx_y_psth/10)];
LW = 1;
subplot(4,1.1,1+1.1)
plot(t_bin1,PSTH{1},'k','linewidth',LW)
axis xy
ylabel('Spike rate [sps]')
set(gca,'ytick',0:20:YL(2));
xlim([0,t_itv1])
ylim(YL)

subplot(4,1.1,1+3*1.1)
plot(t_bin2,PSTH{2},'k','linewidth',LW)
axis xy
xlabel('Time [sec]')
ylabel('Spike rate [sps]')
set(gca,'ytick',0:20:YL(2));
xlim([0,t_itv2])
ylim(YL)

%%
subplot(4,11,11)
barh(1:n_rep,sc{1}/t_itv1,1)
axis ij
ylim(YL)
ylabel('Trials')
set(gca,'ytick',[1,n_rep])
xlabel('Average spike rate [sps]')
ylim(ax_rep)
col = 0.3*[1,1,1];
h = findobj(gca,'Type','patch');
set(h,'FaceColor',col,'EdgeColor',col)


subplot(4,11,3*11)
barh(1:n_rep,sc{2}/t_itv2,1)
axis ij
ylim(YL)
ylabel('Trials')
set(gca,'ytick',[1,n_rep])
xlabel('Average spike rate [sps]')
ylim(ax_rep)
h = findobj(gca,'Type','patch');
set(h,'FaceColor',col,'EdgeColor',col)


subplot(4,11,2*11)
c = 0;
c = c+2;
text(-.2,1-0.1*c,sprintf('- total %d spikes',tot_n_sp{1}))
c = c+2;
text(-.2,1-0.1*c,sprintf('- on average %2.1f [sps]',tot_n_sp{1}/t_itv1/n_rep))
axis off


subplot(4,11,4*11)
c = 0;
c = c+2;
text(-.2,1-0.1*c,sprintf('- total %d spikes',tot_n_sp{2}))
c = c+2;
text(-.2,1-0.1*c,sprintf('- on average %2.1f [sps]',tot_n_sp{2}/t_itv2/n_rep))
axis off

%%

%exportfig(1,'test.eps','height',7,'width',27,'color','rgb')
%orient landscape
%print -dpdf test.pdf
%print('-dpdf',pdf_name)
if t_itv2 == 60
   exportfig(1,pdf_name,'height',7,'width',27,'color','rgb')
elseif t_itv2 == 30
   exportfig(1,pdf_name,'height',7,'width',20,'color','rgb')
end

end
