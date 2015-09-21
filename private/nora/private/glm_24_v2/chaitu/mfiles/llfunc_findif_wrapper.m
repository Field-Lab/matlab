function [f g Hinfo] = llfunc_findif_wrapper(p,pars)

    if(strcmp(pars.basepars.filtermode,'nonsep') && strcmp(pars.basepars.hessmode,'mult'))
        
       [f g Hinfo] = ll_func2_hessmult(p,pars.basepars,pars.stimpars,pars.trainpars); 
    else
        [f g Hinfo] = ll_func2(p,pars.basepars,pars.stimpars,pars.trainpars);
    end