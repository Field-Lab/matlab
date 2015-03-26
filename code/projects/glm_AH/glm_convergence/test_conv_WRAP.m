% Trying Convergence Cut
% AKHEITMAN 2015-02-08

% Scale Convergence to [-0, -1];
% Then choose differe

%%
clear; clc;
BD = NSEM_BaseDirectories;
basefolder   = sprintf('%s/GLM_Convergence_Analysis/fixedSP_rk1_linear_MU_PS_noCP_timekernelCONEMODEL/Analysis_Plots', BD.NSEM_home);

celltypes = [1 2];
exps = [1 2 3 4];


params.st_index = 3; params.st_index_string = '25Percent';
params.taus_unitless = [[linspace(10,1.5,18)],[1./linspace(1,10,19)]];
savedir = sprintf('%s/FitConv_Skip%s', basefolder, params.st_index_string);
if ~exist(savedir,'dir'), mkdir(savedir); end 
i_computation = 1;
i_exp = 1;
i_celltype= 2;
i_cell = 5;
%%
for i_exp = exps
    eval(sprintf('load %s/expTrain_Conv_%d.mat',basefolder,i_exp))
    TC = expTrain_Conv;
    TC.fittedCONV_note = sprintf('Skip First %s of Fit', params.st_index_string);
    exp_nm = TC.exp_nm;
    plotdir = sprintf('%s/%s',savedir,exp_nm);
    if ~exist(plotdir,'dir'), mkdir(plotdir); end
    % LOAD TIMING
    WN.timing_minutes         = TC.WN_fit_minutes;
    NSEM.timing_minutes       = TC.NSEM_fit_minutes;
    %{
    WN.timing.delta_t       =  TC.WN_fit_minutes(end) - TC.WN_fit_minutes(st_index);
    WN.timing.delta_tau     = (TC.WN_fit_minutes(st_index:end) - TC.WN_fit_minutes(st_index)) / WN.timing.delta_t;
    NSEM.timing.delta_t     =  TC.NSEM_fit_minutes(end)- TC.NSEM_fit_minutes(st_index);
    NSEM.timing.delta_tau   = (TC.NSEM_fit_minutes(st_index:end) - TC.NSEM_fit_minutes(st_index)) / NSEM.timing.delta_t;
    %}
    %%
    for i_celltype = celltypes
        % LOAD TRAIN DATA
        if i_celltype == 1 
            cellgroup = TC.ONP;  celltype = 'ONPar';  
            TRAIN = TC.test_train_ONP; XVAL_WN = TC.CONV_WN_ONP; XVAL_NSEM = TC.CONV_NSEM_ONP;  
        end
        if i_celltype == 2
            cellgroup = TC.OFFP; celltype = 'OFFPar'; 
            TRAIN = TC.test_train_OFFP; XVAL_WN = TC.CONV_WN_OFFP; XVAL_NSEM = TC.CONV_NSEM_OFFP;  
        end

        OT = cell(length(cellgroup),1);
        %%
        for i_cell = 1:length(cellgroup)
            cid = cellgroup(i_cell);
            if i_celltype == 1, cell_savename  = sprintf('ONPar_%d',cid);       end
            if i_celltype == 2, cell_savename  = sprintf('OFFPar_%d',cid);      end
            display(sprintf('Working on %s: Cell %s', exp_nm, cell_savename));
            WN.train.raw_vals   = -TRAIN.WN{i_cell}.objval;
            NSEM.train.raw_vals = -TRAIN.NSEM{i_cell}.objval;
            WN.xval.raw_vals    = XVAL_WN{i_cell}.xval_bps;
            NSEM.xval.raw_vals  = XVAL_NSEM{i_cell}.xval_bps;
            
            %%% DEPENDENT ONLY ON WN and NSEM and params
            optimaltime = cell(4,1);
            for i_computation = 1:4
                %%
                % PARAMETERS OF THE COMPUTATION
                tc_tau    = params.taus_unitless; % time constant .. tau (unitless)
                st_index  = params.st_index; %
                clear input
                % GET A LIST OF RAWVALS and FIT_MINUTES 
                if i_computation==1|| i_computation==2, input = WN; end
                if i_computation==3|| i_computation==4, input = NSEM; end
                if i_computation==1|| i_computation==3, rawvals = input.xval.raw_vals(st_index:end);  end
                if i_computation==2|| i_computation==4, rawvals = (input.train.raw_vals(st_index:end))'; end
                delta_minutes = input.timing_minutes(st_index:end) - input.timing_minutes(st_index);
                
                % RESCALE TIME
                delta_tau   = delta_minutes / (delta_minutes(end)); 
                % RESCALE VALUES TO
                unitvals = (rawvals - rawvals(end)) / abs((rawvals(1) - rawvals(end)));
                % DEFINE START OF EXPONENTIAL DECAY FOR VARIOUS TIME CONSTANTS
                startpoint  = 1./(1-exp(-1./tc_tau));
                
                % COMPUTE FIT ERRORS
                l2error = zeros(1,length(tc_tau));
                for i_tau = 1:length(tc_tau)
                    tofit_vals = unitvals - (startpoint(i_tau)-1);
                    fit_vals   = (-startpoint(i_tau)) * exp(-delta_tau/ tc_tau(i_tau));
                    l2error(i_tau) = norm(tofit_vals - fit_vals);
                    %{
                    figure; hold on; 
                    plot(tofit_vals,'r'); plot(fit_vals,'g'); 
                    plot(tofit_vals(1),'k.'); plot(fit_vals(1),'k.');
                    %}
                end
                [min_error,min_ind] = min(l2error) ;
                optimaltime{i_computation}.error         = min_error;
                optimaltime{i_computation}.t_tau         = tc_tau(min_ind);
                optimaltime{i_computation}.t_min         = max(delta_minutes) * tc_tau(min_ind);
                optimaltime{i_computation}.tofit_vals    = unitvals - (startpoint(min_ind)-1);
                optimaltime{i_computation}.fit_vals      =  (-startpoint(min_ind)) * exp(-delta_tau/ tc_tau(min_ind));
                optimaltime{i_computation}.minutes       = delta_minutes;
                optimaltime{i_computation}.skip_minutes  = input.timing_minutes(st_index);
                optimaltime{i_computation}.skip_percent  = input.timing_minutes(st_index) / input.timing_minutes(end);
                
            end
            %%
            % Assign to Cell
            OT{i_cell}.xval_WN_netdecay     = 1/optimaltime{1}.t_tau; 
            OT{i_cell}.xval_NSEM_netdecay   = 1/optimaltime{3}.t_tau;
            OT{i_cell}.train_WN_netdecay    = 1/optimaltime{2}.t_tau; 
            OT{i_cell}.train_NSEM_netdecay  = 1/optimaltime{4}.t_tau; 
            OT{i_cell}.xval_WN_t_min        = optimaltime{1}.t_min; 
            OT{i_cell}.xval_NSEM_t_min      = optimaltime{3}.t_min; 
            OT{i_cell}.train_WN_t_min       = optimaltime{2}.t_min; 
            OT{i_cell}.train_NSEM_t_min     = optimaltime{4}.t_min; 
            OT{i_cell}.skipminutes_WN       = optimaltime{1}.skip_minutes;
            OT{i_cell}.skipminutes_NSEM     = optimaltime{3}.skip_minutes;
            %%
            clf;
            subplot(3,1,1); axis off; hold on;
            c = 0;
            c=c+1; text(-.1, 1-0.1*c,sprintf('Exp: %s  Cell: %s',exp_nm,cell_savename),'interpreter','none');
            c=c+1; text(-.1, 1-0.1*c,sprintf('Column One is Cross Validated Performance (red) in Bits Per Spike')) 
            c=c+1; text(-.1, 1-0.1*c,sprintf('Column Two is Performance on Train Data (green) using the unitless objective function of GLM'))
            c=c+1; text(-.1, 1-0.1*c,sprintf('Skipped First %s of Training Data in Curve Fitting/Analysis', params.st_index_string))
            c=c+1; text(-.1, 1-0.1*c,sprintf('Plot Date %s',datestr(clock)));
            c=c+1; text(-.1, 1-0.1*c,sprintf('Mfile: %s', mfilename('fullpath')),'interpreter','none');
            hold off
            for i_sub = 1:4
                LW_raw = 4; LW_fit = 1.5;
                subplot(3,2,i_sub+2);hold on; set(gca,'fontsize',10);
                xvals       = optimaltime{i_sub}.minutes;
                yvals_tofit = optimaltime{i_sub}.tofit_vals;
                yvals_fit   = optimaltime{i_sub}.fit_vals;
                
                tau   = optimaltime{i_sub}.t_tau;
                t_min = optimaltime{i_sub}.t_min;
                
                switch i_sub
                    case 1
                        heading = 'WN_XVAL'; cstring = 'r';
                    case 2
                        heading = 'WN_TRAIN'; cstring = 'g';
                    case 3
                        heading = 'NSEM_XVAL'; cstring = 'r';
                    case 4
                        heading = 'NSEM_TRAIN'; cstring = 'g';
                end
                title(sprintf('%s: NetDecay: %1.1e, DecayTime_Minutes: %d', heading, 1/tau, round(t_min)),'interpreter','none')
                plot(xvals,yvals_tofit,cstring,'linewidth',LW_raw); plot(xvals,yvals_fit, 'k','linewidth',LW_fit);
                xlim([0 max(xvals)]); xlabel('Fit Minutes'); hold off
            end
            
            orient landscape
            eval(sprintf('print -dpdf %s/%s.pdf',plotdir,cell_savename))
            
        end
        
        if i_celltype == 1, TC.fittedCONV_ONP = OT; end
        if i_celltype == 2, TC.fittedCONV_OFFP = OT; end
    end
    FittedConv = TC;
    eval(sprintf('save %s/fittedCONV_%d.mat FittedConv',savedir,i_exp))
end




%%






%%




