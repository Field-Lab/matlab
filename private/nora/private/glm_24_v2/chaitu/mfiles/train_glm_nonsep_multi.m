function [xstar fstar gstar] = train_glm_nonsep_multi(x0,basepars,stimpars,trainpars,gopts)

gopts = optimset(gopts,'HessMult',@(Hinfo,Y) ll_hess_mult(Hinfo,Y,basepars,stimpars,trainpars),'PrecondBandWidth',1,'MaxPCGIter',basepars.maxCGIter);
if(basepars.maxNonsepIter < optimget(gopts,'MaxIter'))
    gopts = optimset(gopts,'MaxIter',basepars.maxNonsepIter);
end
if (size(stimpars.x,3) == 1)
    pfcn = setup_plotfcns(basepars.filtermode,basepars,stimpars,trainpars);
    gopts = optimset(gopts,'PlotFcns',pfcn);
end

tic; [xstar fstar gstar] = train_glm2_multi(x0,basepars,stimpars,trainpars,gopts,1); toc;
