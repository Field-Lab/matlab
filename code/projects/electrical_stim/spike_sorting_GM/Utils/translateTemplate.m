function templates=translateTemplate(templates,t0,elecs,neuron)

template=templates{neuron};
for e = elecs
    aux                      = template(e,:);
    template(e,1:t0)     = template(e,1);
    template(e,t0+1:end) = aux(1:length(aux)-t0);
end
templates{neuron}=template;