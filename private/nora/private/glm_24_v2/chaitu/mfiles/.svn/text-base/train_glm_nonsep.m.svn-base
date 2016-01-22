function [xstar,fstar,gstar,par_ini,exitflag,output] = train_glm_nonsep(x0,basepars,stimpars,trainpars,gopts,par_ini)
gopts = optimset(...
   gopts,'HessMult',...
   @(Hinfo,Y) ll_hess_mult(Hinfo,Y,basepars,stimpars,trainpars),...
   'PrecondBandWidth',1,...
   'MaxPCGIter',basepars.maxCGIter...
   );

% if(basepars.maxNonsepIter < optimget(gopts,'MaxIter'))
%    gopts = optimset(gopts,'MaxIter',basepars.maxNonsepIter);
% end

%pfcn = setup_plotfcns(basepars.filtermode,basepars,stimpars,trainpars);
%gopts = optimset(gopts,'PlotFcns',pfcn);

if ~exist('par_ini','var')
   par_ini.basepars = basepars;
   par_ini.stimpars = stimpars;
   par_ini.trainpars = trainpars;
   par_ini.gopts = gopts;
end

tic;
[xstar,fstar,gstar,~,exitflag,output] = train_glm2(x0,basepars,stimpars,trainpars,gopts,1);
%[xstar fstar gstar] = train_glm2_lite(x0,basepars,stimpars,trainpars,gopts,1); %edoi
toc;
