%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                         CFD 2 NS Solver                             %%%
%%%                      Author: Senne Hemelaar                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc

%%%========================= Input Variables ===========================%%%
Re_list  = [1000];
N_list   = [15 31 47 55 63];
tol      = 1e-6;

for j = 1:length(Re_list)
for i = 1:length(N_list)
Re = Re_list(j);
N  = N_list(i);
%%%======================== Run the NS Solver ==========================%%%
[results]      = NS_Solver(N, Re, tol);

%%%================= Post Process the Results Vectors ==================%%%
results.xi     = transpose(reshape(results.xi,[N+1,N+1]));
results.ux_xi  = transpose(reshape(results.ux_xi,[N+1,N+1]));
results.uy_xi  = transpose(reshape(results.uy_xi,[N+1,N+1]));
results.tu     = results.Ht11 * results.u;
results        = Reshape_u(results, N);
results        = Reshape_p(results, N);
results        = Reshape_psi(results, N);
results        = Create_Centerline_Results(results);

%%%========================== Save Results =============================%%%
save("Results_N_"+N+"_Re_"+Re+".mat",'results');
end
end









