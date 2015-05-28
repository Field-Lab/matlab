function templates=translateTemplate(templates,t0,elecs,neuron)

template=templates{neuron};
for e = elecs
    aux                      = template(elecs,:);
    template(elecs,1:t0)     = template(elecs,1);
    template(elecs,t0+1:end) = aux(1:length(aux)-t0);
end
templates{neuron}=template;