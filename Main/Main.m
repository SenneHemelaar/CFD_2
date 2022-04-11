%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                         CFD 2 NS Solver                             %%%
%%%                      Author: Senne Hemelaar                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc

%%%========================= Input Variables ===========================%%%
Re       = 1000;
N_list   = [15, 31, 47, 55, 63];
tol      = 1e-6;

for i = 1:length(N_list)
N = N_list(i);
%%%======================== Run the NS Solver ==========================%%%
[results]     = NS_Solver(N, Re, tol);

%%%=================== Reshape the Results Vectors =====================%%%
results.xi    = transpose(reshape(results.xi,[N+1,N+1]));
results.ux_xi = transpose(reshape(results.ux_xi,[N+1,N+1]));
results.uy_xi = transpose(reshape(results.uy_xi,[N+1,N+1]));
results.v     = results.u(length(results.u)/2+1:end);
results.u     = results.u(1:length(results.u)/2);
results.u     = transpose(reshape(results.u,[N+1,N]));
results.v     = transpose(reshape(results.v,[N,N+1]));
results       = Reshape_p(results, N);

%%%========================== Save Results =============================%%%
save("Results_N_"+N+"_Re_"+Re+".mat",'results');
end









