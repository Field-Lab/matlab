function [spks,SU_inp]=generateResp_LNLN_cone_inps(model,movie_c_tf,dt,nTrials)

% calculate subunit activation 
SU_inp = model.cone_to_SU_connection*movie_c_tf;
SU_act = model.fs(SU_inp);

% calculate ganglion cell activation
gang_inp = model.SU_gang_weights'*SU_act;
firing_rate = gang_inp * dt;

% generate spikes
spks = poissrnd(gather(repmat(firing_rate,[nTrials,1])));

spks= double(spks~=0); % make 0 or 1 spike
end