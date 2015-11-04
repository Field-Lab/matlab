function objval = objval_GenLOGISTIC(LOGI_PARAMS,VU, lcif_intoLOGI, lcif_ext, spikebins,t_bin)
MAX     = LOGI_PARAMS(1);
RATE    = LOGI_PARAMS(2);

OFFSET  = log(MAX^VU-1) / RATE;                    
cif_LOGI = MAX ./ ( (1+exp(-RATE*(lcif_intoLOGI- OFFSET))).^(1/VU) );


lcif_LOGI = log(cif_LOGI);
lcif   = lcif_LOGI + lcif_ext;
cif    = exp(lcif);
objval = -( sum( lcif(spikebins) ) - t_bin * sum(cif) );

end