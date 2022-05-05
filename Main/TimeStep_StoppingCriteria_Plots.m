%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%            CFD 2 Time Step and Stopping Criteria Plots              %%%
%%%                      Author: Senne Hemelaar                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; %close all;

%%%========================== Load Results =============================%%%
N        = 31;
Re       = 1000;
dt_list  = [5];
tol_list = logspace(-6,-4,30);

%%%========================= Plot Settings =============================%%%
settings.timestep = 0;
settings.tol      = 0;
settings.RMSD     = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                         Plotting Routines                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (settings.timestep)
t_list = zeros(length(dt_list),1);
for i = 1:length(dt_list)
dt_factor    = dt_list(i);
matfile      = "C:\Users\Senne Hemelaar\OneDrive\Documenten\GitHub\CFD_2\Main\TimeStep_Results\Results_N_"+N+"_Re_"+Re+"_dt_"+dt_factor+"_tol_"+tol+".mat";
results{i}   = load(matfile);
t_list(i)    = results{i}.results.t;
iter_list(i) = results{i}.results.iter;
end

figure(1)
hold on; grid on; box on;
plot(dt_list,t_list,'ks-','Linewidth',1.1)
xlabel('$C_t$','interpreter','latex')
ylabel('CPU Time [s]','interpreter','latex')
f = gcf; filename = "Computing_Time_N_"+N+".png";
file = fullfile("C:\Users\Senne Hemelaar\OneDrive\Documenten\GitHub\CFD_2\Main\Images\TimeStep_Plots", filename);
exportgraphics(f,file,'Resolution',600)

figure(2)
hold on; grid on; box on;
plot(dt_list,iter_list,'ks-','Linewidth',1.1)
xlabel('$C_t$','interpreter','latex')
ylabel('Iterations','interpreter','latex')
f = gcf; filename = "Iterations_N_"+N+".png";
file = fullfile("C:\Users\Senne Hemelaar\OneDrive\Documenten\GitHub\CFD_2\Main\Images\TimeStep_Plots", filename);
exportgraphics(f,file,'Resolution',600)

close all
end

if (settings.tol)
t_list = zeros(length(dt_list),1);
matfile      = "Results_N_"+N+"_Re_"+Re+".mat";
load(matfile);
t_list       = results.t;
iter_list    = results.iter;

% figure(3)
% hold on; grid on; box on;
% plot(results.diff,t_list,'ks-','Linewidth',1.1)
% set(gca, 'XScale', 'log')
% xlim([1e-6 1e-4])
% xlabel('$\epsilon$','interpreter','latex')
% ylabel('CPU Time [s]','interpreter','latex')
% f = gcf; filename = "Tolerance_Computing_Time_N_"+N+".png";
% file = fullfile("C:\Users\Senne Hemelaar\OneDrive\Documenten\GitHub\CFD_2\Main\Images\TimeStep_Plots", filename);
% exportgraphics(f,file,'Resolution',600)

y = linspace(1,results.iter,results.iter-1);

figure(4)
hold on; grid on; box on;
plot(results.diff,y,'k-','Linewidth',1.1)
set(gca, 'XScale', 'log')
xlim([1e-6 1e-4])
xlabel('$\epsilon$','interpreter','latex')
ylabel('Iterations','interpreter','latex')
f = gcf; filename = "Tolerance_Iterations_N_"+N+".png";
file = fullfile("C:\Users\Senne Hemelaar\OneDrive\Documenten\GitHub\CFD_2\Main\Images\TimeStep_Plots", filename);
exportgraphics(f,file,'Resolution',600)
end

if (settings.RMSD)
for i = 1:length(tol_list)
tol = tol_list(i);
dt_factor = 5;
matfile = "C:\Users\Senne Hemelaar\OneDrive\Documenten\GitHub\CFD_2\Main\Tolerance_Results\Results_N_"+N+"_Re_"+Re+"_dt_"+dt_factor+"_tol_"+tol+".mat";
load(matfile);
vec_u{i} = results.vec_u;
vec_u{i} = reshape(vec_u{i},[N^2,1]);
end

for i = 1:length(vec_u)
rmsd(i) = RMSD(vec_u{1}, vec_u{i} ,N);
end

figure(5)
hold on; grid on; box on;
plot(tol_list,rmsd,'k-','Linewidth',1.1)
set(gca, 'XScale', 'log')
xlim([1e-6 1e-4])
xlabel('$\epsilon$','interpreter','latex')
ylabel('RMSE','interpreter','latex')
f = gcf; filename = "Tolerance_Iterations_N_"+N+".png";
file = fullfile("C:\Users\Senne Hemelaar\OneDrive\Documenten\GitHub\CFD_2\Main\Images\TimeStep_Plots", filename);
exportgraphics(f,file,'Resolution',600)
set(gcf,'position',[100,100,750,400])
end



















