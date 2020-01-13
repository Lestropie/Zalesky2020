function mu=calc_mu(N, K, p, v, n1, n2)
  Np=N*p;
  vf=p*v;
  mu=((2*K*v)-v+vf)/(K*(Np+N-n1-n2));
endfunction
