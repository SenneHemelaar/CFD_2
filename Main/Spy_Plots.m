%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                CFD 2 Spy Plots Incidence and Hodge                  %%%
%%%                      Author: Senne Hemelaar                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; %close all;

%%%========================== Load Results =============================%%%
N = 3;
Re = 1000;
matfile = "Results_N_"+N+"_Re_"+Re+".mat";
load(matfile)

%%%========================= Plot Settings =============================%%%
settings.E10   = 0;
settings.tE10  = 0;
settings.E21   = 0;
settings.tE21  = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                         Plotting Routines                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% E10 %%%
if (settings.E10)
E10_pos = results.E10;
E10_neg = results.E10;
for i = 1:24
for j = 1:21
if E10_pos(i,j) < 0.5
    E10_pos(i,j) = 0;
end

if E10_neg(i,j) > -0.5
    E10_neg(i,j) = 0;
end
end
end
figure;
hold on; grid on; box on;
spy(E10_pos,'r')
spy(E10_neg,'b')
set(gca,'xtick',[1:1:21])
set(gca,'ytick',[1:1:24])
xlabel('columns','interpreter','latex')
ylabel('rows','interpreter','latex')
legend('1','-1','interpreter','latex')
set(gcf,'position',[100,100,400,800])
f = gcf; filename = "E10_Spy_Plot.png";
file = fullfile("C:\Users\Senne Hemelaar\OneDrive\Documenten\GitHub\CFD_2\Main\Images\Spy_Plots", filename);
exportgraphics(f,file,'Resolution',600)
end

%%% tE10
if (settings.tE10)
tE10_pos = results.tE10;
tE10_neg = results.tE10;
for i = 1:24
for j = 1:16
if tE10_pos(i,j) < 0.5
    tE10_pos(i,j) = 0;
end

if tE10_neg(i,j) > -0.5
    tE10_neg(i,j) = 0;
end
end
end
figure;
hold on; grid on; box on;
spy(tE10_pos,'r')
spy(tE10_neg,'b')
set(gca,'xtick',[0:1:16])
set(gca,'ytick',[0:1:24])
xlabel('columns','interpreter','latex')
ylabel('rows','interpreter','latex')
legend('1','-1','interpreter','latex')
set(gcf,'position',[100,100,400,550])
f = gcf; filename = "tE10_Spy_Plot.png";
file = fullfile("C:\Users\Senne Hemelaar\OneDrive\Documenten\GitHub\CFD_2\Main\Images\Spy_Plots", filename);
exportgraphics(f,file,'Resolution',600)
end

%%% E21
if (settings.E21)
E21_pos = results.E21;
E21_neg = results.E21;
for i = 1:16
for j = 1:24
if E21_pos(i,j) < 0.5
    E21_pos(i,j) = 0;
end

if E21_neg(i,j) > -0.5
    E21_neg(i,j) = 0;
end
end
end
figure;
hold on; grid on; box on;
spy(E21_pos,'r')
spy(E21_neg,'b')
set(gca,'xtick',[1:1:24])
set(gca,'ytick',[1:1:16])
xlabel('columns','interpreter','latex')
ylabel('rows','interpreter','latex')
legend('1','-1','interpreter','latex')
set(gcf,'position',[100,100,750,400])
f = gcf; filename = "E21_Spy_Plot.png";
file = fullfile("C:\Users\Senne Hemelaar\OneDrive\Documenten\GitHub\CFD_2\Main\Images\Spy_Plots", filename);
exportgraphics(f,file,'Resolution',600)
end

%%% tE21
if (settings.tE21)
tE21_pos = results.tE21;
tE21_neg = results.tE21;
for i = 1:21
for j = 1:24
if tE21_pos(i,j) < 0.5
    tE21_pos(i,j) = 0;
end

if tE21_neg(i,j) > -0.5
    tE21_neg(i,j) = 0;
end
end
end
figure;
hold on; grid on; box on;
spy(tE21_pos,'r')
spy(tE21_neg,'b')
set(gca,'xtick',[1:1:24])
set(gca,'ytick',[1:1:21])
xlabel('columns','interpreter','latex')
ylabel('rows','interpreter','latex')
legend('1','-1','interpreter','latex')
set(gcf,'position',[100,100,900,500])
f = gcf; filename = "tE21_Spy_Plot.png";
file = fullfile("C:\Users\Senne Hemelaar\OneDrive\Documenten\GitHub\CFD_2\Main\Images\Spy_Plots", filename);
exportgraphics(f,file,'Resolution',600)
end

close all;

%%% Ht11











