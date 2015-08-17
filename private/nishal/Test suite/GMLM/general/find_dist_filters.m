function dist = find_dist_filters(filter1,filter2)

[q1,r] = qr(filter1,0);
[q2,r] = qr(filter2,0);

innerprod = q1'*q2;
s = svd(innerprod);
dist = max(sin(acos(s)));

end