function [] = polar_theta_off(h)

rho_labels = {'0' '30' '60' '90' '120' '150' '180' '210' '240' '270' '300' '330'};
for r=1:length(rho_labels)
    delete(findall(h, 'string', rho_labels{r}))
end

end