function ftilde_z_ykp1 = ftilde(K_xp1,B_xp1,K_yp1,B_yp1,gradK,gradB,lambda,fy)


gradterm = sum(sum(gradK.*(K_xp1-K_yp1))) + sum(sum(gradB.*(B_xp1-B_yp1)));
quadterm = norm((K_xp1-K_yp1),'fro') + norm(B_xp1-B_yp1,'fro'); quadterm =quadterm / (2*lambda);
ftilde_z_ykp1 = fy + gradterm + quadterm;
end