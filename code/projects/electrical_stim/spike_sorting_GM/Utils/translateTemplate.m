function templates=translateTemplate(templates,t0,elecs,neuron)
% translateTemplate() shifts the templates to the right  certain
% neuron and electrodes
% 
% inputs:   templates, t0: time units to shift, elecs: electrode indexes in which
% the shift will be done. neuron = neuronindex for the neuron whose
% template will be shifted
% output:  shifted templates
%
% Gonzalo Mena 6/2015 

template=templates{neuron};
for e = elecs
    aux                      = template(e,:);
    template(e,1:t0)     = template(e,1);
    template(e,t0+1:end) = aux(1:length(aux)-t0);
end
templates{neuron}=template;