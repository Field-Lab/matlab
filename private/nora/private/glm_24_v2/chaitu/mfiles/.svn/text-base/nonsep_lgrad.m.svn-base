function lgrad = nonsep_lgrad(basepars,stimpars,trainpars)

lgrad = cell(basepars.Nneurons,1); %zeros(basepars.n*basepars.Mk,size(stimpars.x,2),basepars.Nneurons);
for j=1:basepars.Nneurons
    lgrad{j} = stimpars.dt .* linearfilt_grad_nonsep(0,basepars,stimpars,trainpars,j);
end
