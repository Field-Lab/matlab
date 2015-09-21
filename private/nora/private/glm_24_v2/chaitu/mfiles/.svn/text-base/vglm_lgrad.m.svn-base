function lgrad = vglm_lgrad(basepars,stimpars,trainpars)
% based on nonsep_lgrad.m

% lgrad = cell(basepars.Nneurons,1); %zeros(basepars.n*basepars.Mk,size(stimpars.x,2),basepars.Nneurons);
% for j=1:basepars.Nneurons
%     lgrad{j} = stimpars.dt .* linearfilt_grad_nonsep(0,basepars,stimpars,trainpars,j);
% end

%lgrad = cell(basepars.Nneurons,1); %zeros(basepars.n*basepars.Mk,size(stimpars.x,2),basepars.Nneurons);
lgrad = stimpars.dt .* linearfilt_grad_vglm(basepars,stimpars);
