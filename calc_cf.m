function cf=calc_cf(N, K, p, v, n1, n2)
  Np=N*p;
  vf=p*v;
  mu=calc_mu(N, K, p, v, n1, n2);
  cf=((Np-n2)*mu-vf)^2 + (K-1)*((Np-n2)*mu-v)^2 + K*((N-n1)*mu-v)^2;
endfunction
