function p0 = set_sta_init(rstas,basepars,npars,N,Neff)

p0 = zeros(N*npars,1);

for j=1:N
    [U S V] = svd(reshape(rstas(:,j),basepars.n,basepars.Mk));     % Use SVD to get rank-2 approximation of the STA
    scl = abs(S(1,1));
    offset = (j-1)*npars+1;

    % Spatial filters
    p0(offset+1:offset+basepars.n) = S(1,1)*U(:,1)/scl;
    p0(offset+basepars.n+basepars.Mk+1:offset+2*basepars.n+basepars.Mk) = S(2,2)*U(:,2)/scl;
    % Temporal filters
    p0(offset+basepars.n+1:offset+basepars.n+basepars.Mk) = V(:,1);
    p0(offset+2*basepars.n+basepars.Mk+1:offset+2*(basepars.n+basepars.Mk)) = V(:,2);
    % Postspike filters
    p0(offset+2*(basepars.n+basepars.Mk)+1:offset+2*(basepars.n+basepars.Mk)+basepars.nofilters_postspike) = 0.01 .* randn(basepars.nofilters_postspike,1);
    % Coupling filters
    offset = offset + 2*(basepars.n+basepars.Mk)+basepars.nofilters_postspike;
    p0(offset+1:offset+(Neff-1)*basepars.nofilters_coupling) = 0.01.*randn((Neff-1)*basepars.nofilters_coupling,1);
    
    if 0 % Visualize Rank-2 approximation
        [s1 t1 s2 t2] = get_sep_filters(p0((j-1)*npars+1:j*npars),basepars.n,basepars.Mk);
        sta_r2 = s1*t1' + s2*t2'; % Rank-2 approximation of the sta
        sta_r2 = sta_r2 - min(sta_r2(:));
        sta_r2 = sta_r2./max(sta_r2(:));
        play_sta(sta_r2,sqrt(basepars.n),sqrt(basepars.n),basepars.Mk);
    end
end