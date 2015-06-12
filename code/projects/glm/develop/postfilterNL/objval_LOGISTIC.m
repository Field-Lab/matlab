function objval = objval_LOGISTIC(LOGI_PARAMS, lcif_intoLOGI, lcif_ext, spikebins,t_bin)
MAX     = LOGI_PARAMS(1);
RATE    = LOGI_PARAMS(2);


OFFSET  = log(MAX-1) / RATE;                    

lcif_LOGI = log(MAX ./ (1 + exp(-RATE * (lcif_intoLOGI- OFFSET) )));

lcif   = lcif_LOGI + lcif_ext;
cif    = exp(lcif);
objval = -( sum( lcif(spikebins) ) - t_bin * sum(cif) );
end