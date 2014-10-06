clear; clear all; close all; clc

 
cellselectiontype = 'shortlist';
exptests = [1 2 3 4];
modulation  = 'model_parameters';

%{
mod_version_base = 'ONOFFvsLin_8pix_1e4_8pix';
purpose_base = 'Examine effect of ON/OFF Hard rectification channels, WITH cone model, fixed spatial filter STA';
comp_params{1}.name = 'single linear channel';
comp_params{2}.name = 'hard rect ON and Off channel';
comp_params{1}.filtertype = 'fixedSP';
comp_params{2}.filtertype = 'OnOff_hardrect_fixedSP_STA';

mod_version_base = 'ONOFFvsLin_8pix_ID_8pix';
purpose_base = 'Examine effect of ON/OFF Hard rectification channels, withOUT cone model, fixed spatial filter STA';
comp_params{1}.name = 'single linear channel';
comp_params{2}.name = 'hard rect ON and Off channel';
comp_params{1}.filtertype = 'fixedSP';
comp_params{2}.filtertype = 'OnOff_hardrect_fixedSP_STA';
%}

mod_version_base = 'fixedSPvsrk2_8pix_ID_8pix';
purpose_base = 'Examine effect of rank-2, withOUT cone model, fixedSP uses STA for spatial filter';
comp_params{1}.name = 'fixedSP';
comp_params{2}.name = 'rk2';
comp_params{1}.filtertype = 'fixedSP';
comp_params{2}.filtertype = 'rk2';

mod_version_base = 'fixedSPvsrk2_8pix_1e4_8pix';
purpose_base = 'Examine effect of rank-2, with cone model, fixedSP uses STA for spatial filter';
comp_params{1}.name = 'fixedSP';
comp_params{2}.name = 'rk2';
comp_params{1}.filtertype = 'fixedSP';
comp_params{2}.filtertype = 'rk2';

mod_version_base = 'conemodel_8pix_1e4_8pix';
purpose_base = 'Examine effect of  cone model, fixed spatial filter STA';
comp_params{1}.name = 'no cone model';
comp_params{2}.name = 'with cone model';
comp_params{1}.filtertype = 'fixedSP';
comp_params{2}.filtertype = 'fixedSP';



for i_test= 1:2
    clear fullmodelcomparison
    if i_test == 1
        eval(sprintf('load WN_%s/shortlist_crossprep.mat', mod_version_base))
        mcWN   = fullmodelcomparison;
    end
    if i_test == 2
        eval(sprintf('load NSEM_%s/shortlist_crossprep.mat', mod_version_base))
        mcNSEM = fullmodelcomparison;
    end
end
clear fullmodelcomparison


%%
clf
plotmfile = mfilename('fullpath');
subplot(6,1,4);
axis off
set(gca, 'fontsize', 12)
c = 0;
c = c+1;
text(-.1, .8-0.1*c,sprintf('Purpose is to: %s', purpose_base ), 'FontSize',8);
c = c+1;
text(-.1, .8-0.1*c,sprintf('UOP: Unconditioned Optimal Prediction, ie. the mean rate of the raster'), 'FontSize',8);
c = c+1;
text(-.1, .8-0.1*c,sprintf('"Bits" refers to the (LogProb of Model - LogProb of Uniform Firing Rate)'), 'FontSize',8);
c = c+1;
text(-.1, .8-0.1*c,sprintf('Symbol * is for OFF Parasol, Symbol . is for ON Parasol'), 'FontSize',8);
c = c+1;
text(-.1, .8-0.1*c,sprintf('If bits are negative, consider model fail , replotted as zero'),'FontSize',8);
c = c+1;
text(-.1, .8-0.1*c,sprintf('Plot Date %s',datestr(clock)),'FontSize',8);
c= c+1;
text(-.1, .8-0.1*c,sprintf('Plot mfile: %s',plotmfile ),'FontSize',8 );


MS = 16;

for i_stim = 1:2
    amax= [];
    amin = [];
    for i_exp = exptests 
        
        

        if i_stim ==1 
            compare = mcWN{i_exp};
            subplot(6,2,[3 5]); hold on;
            xlabel(compare.comparenames{1}); ylabel(compare.comparenames{2});
            set(gca, 'fontsize', 10);
            title('WN GLMperformance/UOPperformance');     
        end
        if i_stim ==2
            compare = mcNSEM{i_exp}; 
            subplot(6,2,[4 6]); hold on;
            xlabel(compare.comparenames{1}); ylabel(compare.comparenames{2});
            set(gca, 'fontsize', 10);
            title('NSEM  GLMperformance/UOPperformance');  
               
        end
        
        measure =  (compare.logprob_glm_bpspike ./ compare.logprob_uop_bpspike);
        
        ONP =  find(compare.ONPar);
        OFFP = find(compare.OFFPar);

        modelfails = find(measure<=0);
        measure(modelfails) = 0 ;
        amax = [amax max(measure(:,:))];            
        newmin = min(min(min(measure(:,ONP))) , min(min(measure(:,OFFP))));
        amin = [amin newmin];

        switch i_exp
            case 1
                plot(  measure(1,OFFP),  measure(2,OFFP)  ,'r*', 'markersize' , MS ); 
                plot(  measure(1,ONP),   measure(2,ONP)  ,'r.', 'markersize' , MS ); 
            case 2
                plot(  measure(1,OFFP),  measure(2,OFFP)  ,'g*', 'markersize' , MS ); 
                plot(  measure(1,ONP),  measure(2,ONP)  ,'g.', 'markersize' , MS ); 
            case 3
                plot(  measure(1,OFFP),  measure(2,OFFP)  ,'b*', 'markersize' , MS ); 
                plot(  measure(1,ONP),  measure(2,ONP)  ,'b.', 'markersize' , MS ); 
            case 4
                plot(  measure(1,OFFP),  measure(2,OFFP)  ,'k*', 'markersize' , MS ); 
                plot(  measure(1,ONP),  measure(2,ONP)  ,'k.', 'markersize' , MS ); 
        end
            
    end
	zmax = max(amax);
	zmin = 0% min(amin);
	plot(linspace(zmin,zmax,100), linspace(zmin,zmax,100),'k');
	xlim([zmin zmax]) ; ylim([zmin zmax]);
end

subplot(6,1,1);
axis off
set(gca, 'fontsize', 12)
c = 0;


text(-.1, .8,sprintf('%s',purpose_base),'FontSize',10 );
text(-.1, .6, sprintf('Upper Plots Compare Fits qualities for WN (left) and NSEM(right)'),'FontSize',10);
text(-.1, .4, sprintf('Lower Plots Compare Change in WN to NSEM quality ratio'),'FontSize',10);



zmin = 0;
zmax = 0;
for i_fit = 1:2
    amax= [];
    amin = [];
    for i_exp = exptests 

        if i_fit == 1 
            subplot(6,2,[9 11]); hold on; 
            xlabel('WN fit and test'); ylabel('NSEM fit and test');
            set(gca, 'fontsize', 10);
            title(sprintf('WN vs NSEM: %s',comp_params{1}.name));
        end
        if i_fit == 2
            compare = mcNSEM{i_exp}; 
            subplot(6,2,[10 12]); hold on;
            xlabel('WN fit and test'); ylabel('NSEM fit and test');
            set(gca, 'fontsize', 10);
            title(sprintf('WN vs NSEM: %s',comp_params{2}.name)); 
        end
       
        measureWN   = (mcWN{i_exp}.logprob_glm_bpspike(i_fit,:) ./ mcWN{i_exp}.logprob_uop_bpspike(i_fit,:));
        measureNSEM   = (mcNSEM{i_exp}.logprob_glm_bpspike(i_fit,:) ./ mcNSEM{i_exp}.logprob_uop_bpspike(i_fit,:));
        
        measure = [measureWN ; measureNSEM];
        
        ONP =  find(mcWN{i_exp}.ONPar);
        OFFP = find(mcWN{i_exp}.OFFPar);

        modelfails = find(measure<=0);
        measure(modelfails) = 0 ;
        amax = [amax max(measure(:,:))];            
        newmin = min(min(min(measure(:,ONP))) , min(min(measure(:,OFFP))));
        amin = [amin newmin];

       switch i_exp
            case 1
                plot(  measure(1,OFFP),  measure(2,OFFP)  ,'r*', 'markersize' , MS ); 
                plot(  measure(1,ONP),   measure(2,ONP)  ,'r.', 'markersize' , MS ); 
            case 2
                plot(  measure(1,OFFP),  measure(2,OFFP)  ,'g*', 'markersize' , MS ); 
                plot(  measure(1,ONP),  measure(2,ONP)  ,'g.', 'markersize' , MS ); 
            case 3
                plot(  measure(1,OFFP),  measure(2,OFFP)  ,'b*', 'markersize' , MS ); 
                plot(  measure(1,ONP),  measure(2,ONP)  ,'b.', 'markersize' , MS ); 
            case 4
                plot(  measure(1,OFFP),  measure(2,OFFP)  ,'k*', 'markersize' , MS ); 
                plot(  measure(1,ONP),  measure(2,ONP)  ,'k.', 'markersize' , MS ); 
        end
            
    end
	zmax = max(max(amax),zmax);
end

for i_fit = 1:2
    if i_fit == 1, subplot(6,2,[9 11]);  end
    if i_fit == 2, subplot(6,2,[10 12]); end
	plot(linspace(zmin,zmax,100), linspace(zmin,zmax,100),'k');
	xlim([zmin zmax]) ; ylim([zmin zmax]);
end


%
orient tall
eval(sprintf('print -dpdf %s_%s.pdf',cellselectiontype,mod_version_base));









