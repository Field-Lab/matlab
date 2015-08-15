function NLOutput = LogFixMu_withPS_fmincon(NL_Input,home_spbins,t_bin)
% 2015-06-24  AKHeitman

% NL_Input needs: y_int (value with no drive)
%                input_test(corresponding input to the NL for fit data)
%                input_fit(corresponding input to the NL for test data) 

optim_struct = optimset(...
                    'derivativecheck','off','diagnostics','off',...  % 
                    'display','iter','funvalcheck','off',... 
                    'MaxIter',100,'TolFun',10^(-6),'TolX',10^(-9) );
lowerbound = [NL_Input.y_int+1 .1];
upperbound = [1000         100];
LOGI_Params0 = [100,1];

if ~isfield(NL_Input,'fit_lcif_additionaldrive')
    [LOGI_Params_Opt, new_objval, eflag, output] = fmincon(@(LOGI_Params) subR_objval_LOGISTIC...
        (LOGI_Params, NL_Input.y_int,NL_Input.input_fit, 0, home_spbins,t_bin),...
        LOGI_Params0,[],[],[],[],lowerbound,upperbound,[],optim_struct);
    
    [~, lcif_LOGI_crossvaltest] =  subR_objval_LOGISTIC(LOGI_Params_Opt,...
    NL_Input.y_int,NL_Input.input_test, 0, [],t_bin); 
    
    [~, lcif_LOGI_fit] =  subR_objval_LOGISTIC(LOGI_Params_Opt,...
    NL_Input.y_int,NL_Input.input_fit, 0, home_spbins,t_bin); 

else
    [LOGI_Params_Opt, new_objval, eflag, output] = fmincon(@(LOGI_Params) subR_objval_LOGISTIC...
        (LOGI_Params, NL_Input.y_int,NL_Input.input_fit, NL_Input.fit_lcif_additionaldrive, home_spbins,t_bin),...
        LOGI_Params0,[],[],[],[],lowerbound,upperbound,[],optim_struct);
    
    [~, lcif_LOGI_crossvaltest] =  subR_objval_LOGISTIC(LOGI_Params_Opt,...
    NL_Input.y_int,NL_Input.input_test, 0, [],t_bin); 
    
    [~, lcif_LOGI_fit] =  subR_objval_LOGISTIC(LOGI_Params_Opt,...
    NL_Input.y_int,NL_Input.input_fit, NL_Input.fit_lcif_additionaldrive, home_spbins,t_bin);  
end






NLOutput.note_metric = 'Logistic with Null Rate fixed to tonic drive. Driven by linear filter output,normalized to std 1';
NLOutput.param_string = sprintf('2 Params fit with fmincon,  Optimal Max Rate: %1.1e, Optimal Slope: %1.1e',...
    LOGI_Params_Opt(1), LOGI_Params_Opt(2) );
NLOutput.maxrate = LOGI_Params_Opt(1);
NLOutput.slope   = LOGI_Params_Opt(2);
NLOutput.fmincon_output = output;
NLOutput.eflag = eflag;
NLOutput.crossvaltest_finalrate = exp(lcif_LOGI_crossvaltest);
NLOutput.new_objval = new_objval;
NLOutput.fit_rate = exp(lcif_LOGI_fit);
NLOutput.codename = mfilename('fullpath');
clear lcif_LOGI_crossvaltest dummy lowerbound upperbound LOGI_Params0
% Check Code
%{
lcif   = baseGLM.lcif_fit.stim + baseGLM.lcif_fit.mu; 
cif    = exp(lcif);
objval = -( sum( lcif(home_spbins) ) - t_bin * sum(cif) );
pars = [10,1];
objval = objval_LOGISTIC(pars,NL_Input.y_int,NL_Input.input_fit, 0, home_spbins,t_bin)
%}
end

function [objval lcif_LOGI] = subR_objval_LOGISTIC(LOGI_PARAMS, Y_INT, lcif_intoLOGI, lcif_ext, spikebins,t_bin)
MAX     = LOGI_PARAMS(1);
RATE    = LOGI_PARAMS(2);


OFFSET  = log( (MAX/Y_INT) - 1  ) / RATE;                    

lcif_LOGI = log(MAX ./ (1 + exp(-RATE * (lcif_intoLOGI- OFFSET) )));

lcif   = lcif_LOGI + lcif_ext;
cif    = exp(lcif);
objval = -( sum( lcif(spikebins) ) - t_bin * sum(cif) );
end