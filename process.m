% Reproduce a number of steps that were presented in:
%
% Zalesky A, Sarwar T, Ramamohanarao K
% A cautionary note on the use of SIFT in pathological connectomes
% Magn Reson Med. 2020 Mar;83(2):791-794. doi: 10.1002/mrm.28037
%
% , but presenting additional data derived from those experiments.

% Additional files required:
%   calc_cf.m
%   calc_mu.m
%   calc_n1.m

% Note: Figure positioning / font sizes / plot widths are optimised for a 4K display



clear all;
close all;

% constants throughout experiment
v=1.0;
n2=0;



% draw raw streamline counts as points rather than lines, since
%   they are not derivatives of the experiment
N=40;
N1=repmat(N-n2, 1, N);
N2=1:N;
N1_perc=-100.0*(N1-N)/N;
N2_perc=-100.0*(N2-N)/N;

figure;
subplot(2, 9, 1);
hold on
title("N ; N x p \n(streamlines counts)", "FontSize", 20)
axis("Tic", "Labely")
set(gca, "XTickMode", "manual")
set(gca, "XTick", [1 10 100])
set(gca, "FontSize", 28)
ylabel("Percentage Change", "FontSize", 40)
axis([1 100 0 100])
semilogx(1:N, N1_perc, "color", "blue", "marker", "+", "markersize", 15, "linestyle", "none")
semilogx(1:N, N2_perc, "color", "red", "marker", "+", "markersize", 15, "linestyle", "none")

subplot(2, 9, 10);
hold on
set(gca, "XTickMode", "manual")
set(gca, "XTick", [1 10 100])
xticklabels([{"10^0"} {"10^1"} {"10^2"}])
set(gca, "FontSize", 28)
ylabel("Connectivity", "FontSize", 40)
axis([1 100 0 40])
semilogx(1:N, N1, "color", "blue", "marker", "+", "markersize", 15, "linestyle", "none")
semilogx(1:N, N2, "color", "red", "marker", "+", "markersize", 15, "linestyle", "none")



column_index=2;
for N=[40,1000];

  Nrange=1:N;

  % generate the appropriate x-axis ticks based on the number of streamlines
  xaxis_max_log=ceil(log10(N));
  xaxis_max=10^xaxis_max_log;
  xaxis_tick_powers=0:1:xaxis_max_log;
  xaxis_tick_values=10.^xaxis_tick_powers;
  xaxis_tick_strings=repmat({""}, 1, length(xaxis_tick_powers));
  for i=1:1:length(xaxis_tick_strings)
    xaxis_tick_strings(i)=strcat("10^",num2str(xaxis_tick_powers(i)));
  end

  for K=[5,10,40];
    c1_exact=zeros(1,N);
    c2_exact=zeros(1,N);
    c1_round=zeros(1,N);
    c2_round=zeros(1,N);
    c1_lower=zeros(1,N);
    c2_lower=zeros(1,N);
    c1_upper=zeros(1,N);
    c2_upper=zeros(1,N);
    for Np=Nrange;
      p=Np/N;
      n1=calc_n1(N, K, p, v);
      mu_exact=calc_mu(N, K, p, v, n1, n2);
      c1_exact(1,Np)=(N-n1)*mu_exact;
      c2_exact(1,Np)=Np*mu_exact;

      % find range of floating-point solutions that would have mapped
      %   onto the same integer solution
      n1_floor=floor(n1);
      n1_ceil=ceil(n1);
      cf_floor=calc_cf(N, K, p, v, n1_floor, n2);
      cf_ceil=calc_cf(N, K, p, v, n1_ceil, n2);
      if cf_ceil < cf_floor
        n1_round=n1_ceil;
        n1_lower=n1_floor;
        if calc_cf(N, K, p, v, n1_ceil+1, n2) < cf_ceil
          n1_upper=n1_ceil;
        else
          n1_upper=n1_ceil+1;
        endif
      else
        n1_round=n1_floor;
        if calc_cf(N, K, p, v, n1_floor-1, n2) < cf_floor
          n1_lower=n1_floor;
        else
          n1_lower=n1_floor-1;
        endif
        n1_upper=n1_ceil;
      endif

      mu_round=calc_mu(N, K, p, v, n1_round, n2);
      c1_round(1,Np)=(N-n1_round)*mu_round;
      c2_round(1,Np)=Np*mu_round;
      mu_lower=calc_mu(N, K, p, v, n1_lower, n2);
      c1_lower(1,Np)=(N-n1_lower)*mu_lower;
      c2_lower(1,Np)=Np*mu_lower;
      mu_upper=calc_mu(N, K, p, v, n1_upper, n2);
      c1_upper(1,Np)=(N-n1_upper)*mu_upper;
      c2_upper(1,Np)=Np*mu_upper;
    endfor

    c1_exact_perc=-100.0*(c1_exact-v)/v;
    c2_exact_perc=-100.0*(c2_exact-v)/v;
    c1_round_perc=-100.0*(c1_round-v)/v;
    c2_round_perc=-100.0*(c2_round-v)/v;
    c1_lower_perc=-100.0*(c1_lower-v)/v;
    c2_lower_perc=-100.0*(c2_lower-v)/v;
    c1_upper_perc=-100.0*(c1_upper-v)/v;
    c2_upper_perc=-100.0*(c2_upper-v)/v;

    subplot(2, 9, column_index);

    title(strcat("N=", num2str(N), ";\nK=", num2str(K), " (", num2str(2*K), "mm)"), "FontSize", 20)
    axis("Tic", "Labely")
    set(gca, "XTickMode", "manual")
    set(gca, "XTick", xaxis_tick_values)
    set(gca, "FontSize", 28)

    % specify display range along y-axis to match Zalesky et al., 2019
    % for data omitted from Zalesky et al., 2019, use the same range as the
    %   corresponding sub-figure with a smaller streamline count
    if N==40
      if K==5
        yrange=[-10; 30];
      elseif K==10
        yrange=[-5; 15];
      elseif K==40
        yrange=[0; 1.5];
      endif
    elseif N==1000
      if K==5
        yrange=[-10; 30];
      elseif K==10
        yrange=[-5; 15];
      elseif K==40
        yrange=[-2; 4];
      endif
    endif
    axis(cat(1, [1; xaxis_max], yrange))

    hold on

    patch([Nrange, Nrange(end:-1:1)]', [c1_upper_perc, c1_lower_perc(end:-1:1)]', "blue", "linestyle", ":", "edgecolor", "blue", "facealpha", 0.1)
    patch([Nrange, Nrange(end:-1:1)]', [c2_upper_perc, c2_lower_perc(end:-1:1)]', "red", "linestyle", ":", "edgecolor", "red", "facealpha", 0.1)

    semilogx(Nrange, c1_round_perc, "color", "blue", "linewidth", 6)
    semilogx(Nrange, c2_round_perc, "color", "red", "linewidth", 6)
    semilogx(Nrange, c1_exact_perc, "color", "blue", "linestyle", "--", "linewidth", 4)
    semilogx(Nrange, c2_exact_perc, "color", "red", "linestyle", "--", "linewidth", 4)



    subplot(2, 9, 9+column_index)
    set(gca, "XTickMode", "manual")
    set(gca, "XTick", xaxis_tick_values)
    set(gca, "FontSize", 28)

    if N==40
      if K==5
        yrange=[0.7; 1.1];
      elseif K==10
        yrange=[0.85; 1.05];
      elseif K==40
        yrange=[0.985; 1.00];
      else
        error
      endif
    elseif N==1000
      if K==5
        yrange=[0.7; 1.1];
      elseif K==10
        yrange=[0.85; 1.05];
      elseif K==40
        yrange=[0.96; 1.02];
      else
        error
      endif
    else
      error
    endif
    axis(cat(1, [1; xaxis_max], yrange))
    xticklabels(xaxis_tick_strings)

    hold on

    patch([Nrange, Nrange(end:-1:1)]', [c1_upper, c1_lower(end:-1:1)]', "blue", "linestyle", ":", "edgecolor", "blue", "facealpha", 0.1)
    patch([Nrange, Nrange(end:-1:1)]', [c2_upper, c2_lower(end:-1:1)]', "red", "linestyle", ":", "edgecolor", "red", "facealpha", 0.1)

    semilogx(Nrange, c1_round, "color", "blue", "linewidth", 6)
    semilogx(Nrange, c2_round, "color", "red", "linewidth", 6)
    semilogx(Nrange, c1_exact, "color", "blue", "linestyle", "--", "linewidth", 4)
    semilogx(Nrange, c2_exact, "color", "red", "linestyle", "--", "linewidth", 4)

    column_index = column_index + 1;
  endfor
endfor






N=40;
Nrange=1:N;
Klist=[5 10 20 40];
colors=[[0.59 0.00 0.20]; [0.87 0.22 0.32]; [0.94 0.40 0.22]; [0.99 0.68 0.37]];
c1_exact_perc=zeros(4,N);
c1_exact=zeros(4,N);
c2_exact_perc=zeros(4,N);
c2_exact=zeros(4,N);
for Kindex=1:4;
  K=Klist(Kindex);
  for Np=1:N;
    p=Np/N;
    n1=calc_n1(N, K, p, v);
    mu_exact=calc_mu(N, K, p, v, n1, n2);
    c1_exact(Kindex,Np)=(N-n1)*mu_exact;
    c2_exact(Kindex,Np)=Np*mu_exact;
  endfor
endfor

c1_exact_perc=-100.0*(c1_exact-v)/v;
c2_exact_perc=-100.0*(c2_exact-v)/v;

subplot(2, 9, 8);
title("Non-Integer\nSolution")
axis("Tic", "Labely")
set(gca, "XTickMode", "manual")
set(gca, "XTick", xaxis_tick_values)
set(gca, "FontSize", 28)
axis([1 100 0 20])
hold on
for Kindex=1:4;
  semilogx(Nrange, c2_exact_perc(Kindex,:), "color", colors(Kindex,:), "linewidth", 6, "linestyle", "--")
  semilogx(Nrange, c1_exact_perc(Kindex,:), "color", [0.54 0.63 0.80], "linewidth", 6, "linestyle", "--")
endfor

subplot(2, 9, 17);
set(gca, "XTickMode", "manual")
set(gca, "XTick", [1 10 100])
xticklabels([{"10^0"} {"10^1"} {"10^2"}])
set(gca, "FontSize", 28)
axis([1 100 0.75*v v])
hold on
for Kindex=1:4;
  semilogx(Nrange, c2_exact(Kindex,:), "color", colors(Kindex,:), "linewidth", 6, "linestyle", "--")
  semilogx(Nrange, c1_exact(Kindex,:), "color", [0.54 0.63 0.80], "linewidth", 6, "linestyle", "--")
endfor



FAref=0.8;
c1_FA_perc=zeros(4,N);
c1_FA=zeros(4,N);
c2_FA_perc=zeros(4,N);
c2_FA=zeros(4,N);
for Kindex=1:4;
  K=Klist(Kindex);
  for Np=1:N;
    p=Np/N;
    c1_FA(Kindex,Np)=(K*FAref)/K;
    c2_FA(Kindex,Np)=(((K-1)*FAref)+(FAref*p))/K;
  endfor
endfor

c1_FA_perc=-100.0*(c1_FA-FAref)/FAref;
c2_FA_perc=-100.0*(c2_FA-FAref)/FAref;

subplot(2, 9, 9);
hold on
title("Mean FA")
axis("Tic", "Labely")
set(gca, "XTickMode", "manual")
set(gca, "XTick", [1 10 100])
set(gca, "FontSize", 28)
axis([1 100 0 20])
for Kindex=1:4;
  semilogx(Nrange, c2_FA_perc(Kindex,:), "color", colors(Kindex,:), "linewidth", 6)
  semilogx(Nrange, c1_FA_perc(Kindex,:), "color", [0.54 0.63 0.80], "linewidth", 6)
endfor

subplot(2, 9, 18);
hold on
set(gca, "XTickMode", "manual")
set(gca, "XTick", [1 10 100])
xticklabels([{"10^0"} {"10^1"} {"10^2"}])
set(gca, "FontSize", 28)
axis([1 100 0.75*FAref FAref])
for Kindex=1:4;
  semilogx(Nrange, c2_FA(Kindex,:), "color", colors(Kindex,:), "linewidth", 6)
  semilogx(Nrange, c1_FA(Kindex,:), "color", [0.54 0.63 0.80], "linewidth", 6)
endfor






ha=get(gcf,'children');
set(ha(17),'position',[0.03 0.04 0.08 0.40])
set(ha(15),'position',[0.17 0.04 0.08 0.40])
set(ha(13),'position',[0.27 0.04 0.08 0.40])
set(ha(11),'position',[0.37 0.04 0.08 0.40])
set(ha(9),'position',[0.47 0.04 0.08 0.40])
set(ha(7),'position',[0.57 0.04 0.08 0.40])
set(ha(5),'position',[0.67 0.04 0.08 0.40])
set(ha(3),'position',[0.81 0.04 0.08 0.40])
set(ha(1),'position',[0.91 0.04 0.08 0.40])
set(ha(18),'position',[0.03 0.51 0.08 0.40])
set(ha(16),'position',[0.17 0.51 0.08 0.40])
set(ha(14),'position',[0.27 0.51 0.08 0.40])
set(ha(12),'position',[0.37 0.51 0.08 0.40])
set(ha(10),'position',[0.47 0.51 0.08 0.40])
set(ha(8),'position',[0.57 0.51 0.08 0.40])
set(ha(6),'position',[0.67 0.51 0.08 0.40])
set(ha(4),'position',[0.81 0.51 0.08 0.40])
set(ha(2),'position',[0.91 0.51 0.08 0.40])

