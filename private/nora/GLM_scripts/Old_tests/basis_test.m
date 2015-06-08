GLMPars = GLMParams;

bin_size = 0.0083 / GLMPars.bins_per_frame;
basis_params  = GLMPars.spikefilters.cp;
cp_basis      = prep_spikefilterbasisGP(basis_params,bin_size);

basis_params  = GLMPars.spikefilters.ps;
ps_basis      = prep_spikefilterbasisGP(basis_params,bin_size);

figure;
plot(cp_basis)

figure;
plot(ps_basis)