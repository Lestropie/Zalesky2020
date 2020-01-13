function n1=calc_n1(N, K, p, v)
  Np=N*p;
  vf=p*v;
  n1=N-((K*Np*v)/(vf-v+(K*v)));
endfunction
