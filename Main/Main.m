%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                         CFD 2 NS Solver                             %%%
%%%                      Author: Senne Hemelaar                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc

%%%========================= Input Variables ===========================%%%
Re_list   = [1000];
N_list    = [31];
dt_list   = [5];
tol_list  = [1e-6];

%%%========================= Loop Over Inputs ==========================%%%
for m = 1:length(tol_list)
for k = 1:length(dt_list)
for j = 1:length(Re_list)
for i = 1:length(N_list)
tic
Re        = Re_list(j);
N         = N_list(i);
dt_factor = dt_list(k);
tol       = tol_list(m);

%%%======================== Run the NS Solver ==========================%%%
[results]      = NS_Solver(N, Re, tol, dt_factor);

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
if length(dt_list) > 2
    save("C:\Users\Senne Hemelaar\OneDrive\Documenten\GitHub\CFD_2\Main\TimeStep_Results\Results_N_"...
         +N+"_Re_"+Re+"_dt_"+dt_factor+"_tol_"+tol+".mat",'results');
elseif length(tol_list) >2
    save("C:\Users\Senne Hemelaar\OneDrive\Documenten\GitHub\CFD_2\Main\Tolerance_Results\Results_N_"...
         +N+"_Re_"+Re+"_dt_"+dt_factor+"_tol_"+tol+".mat",'results');
else
    save("Results_N_"+N+"_Re_"+Re+".mat",'results')
end

end
end
end
end









