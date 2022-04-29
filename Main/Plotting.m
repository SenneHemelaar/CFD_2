%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                           CFD 2 Plotter                             %%%
%%%                      Author: Senne Hemelaar                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all

%%%========================= Plot Settings =============================%%%
settings.pressure_contour_plot   = 0;
settings.vorticity_contour_plot  = 0;
settings.streamline_contour_plot = 0;

settings.velocity_line_plot_Cx   = 1;
settings.pressure_line_plot_Cx   = 1;
settings.vorticity_line_plot_Cx  = 1;

settings.velocity_line_plot_Cy   = 1;
settings.pressure_line_plot_Cy   = 1;
settings.vorticity_line_plot_Cy  = 1;

%%%========================== Load Results =============================%%%
N_list = [15 31 63];
for i = 1:length(N_list)
N = N_list(i);
Re = 1000;
matfile = "Results_N_"+N+"_Re_"+Re+".mat";
load(matfile)

%%%======================= Load Reference Data =========================%%%
load('reference_data_x05.txt');
ref_Cx.y    = reference_data_x05(:,1);
ref_Cx.uref = reference_data_x05(:,2);
ref_Cx.u    = reference_data_x05(:,3);
ref_Cx.p    = reference_data_x05(:,4);
ref_Cx.xi   = reference_data_x05(:,5);

load('reference_data_y05.txt');
ref_Cy.x    = reference_data_y05(:,1);
ref_Cy.vref = reference_data_y05(:,2);
ref_Cy.v    = reference_data_y05(:,3);
ref_Cy.p    = reference_data_y05(:,4);
ref_Cy.xi   = reference_data_y05(:,5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                         Plotting Routines                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%========================== Contour Plots ============================%%%

%%% PRESSURE %%%
if (settings.pressure_contour_plot)
[x,y] = meshgrid(results.x(2:end-1),results.x(2:end-1));
figure;
[C,h] = contour(x,y,results.p,[-0.002 0.0 0.02 0.05 0.07...
        0.09 0.11 0.12 0.17 0.3],'LineColor','black');
v = [0.0 0.02 0.05 0.09 0.11 0.17 0.3];
clabel(C,h,v);
xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')
f = gcf; filename = "Pressure_Contour_Re_"+Re+"_N_"+N+".png";
file = fullfile("C:\Users\Senne Hemelaar\OneDrive\Documenten\GitHub\CFD_2\Main\Images\Contours", filename);
exportgraphics(f,file,'Resolution',600)
end

%%% VORTICITY %%%
if (settings.vorticity_contour_plot)
[x,y] = meshgrid(results.tx,results.tx);
figure;
[C,h] = contour(x,y,results.xi,[-3.0 -2.0 -1.0 -0.5 0.0...
        0.5 1.0 2.0 3.0 4.0 5.0],'LineColor','black');
v = [-3.0 -2.0 -1.0 0.0 1.0 2.0 3.0 5.0];
clabel(C,h,v)
xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')
f = gcf; filename = "Vorticity_Contour_Re_"+Re+"_N_"+N+".png";
file = fullfile("C:\Users\Senne Hemelaar\OneDrive\Documenten\GitHub\CFD_2\Main\Images\Contours", filename);
exportgraphics(f,file,'Resolution',600)
end

%%% STREAMLINES %%%
if (settings.streamline_contour_plot)
[x,y] = meshgrid(results.tx,results.tx);
figure;
[C,h] = contour(x,y,results.psi,[0.1175 0.115 0.11 0.1 9e-2 7e-2...
        5e-2 3e-2 1e-2 1e-4 1e-5 1e-10 0.0 -1e-6 -1e-5 -5e-5 -1e-4...
        -2.5e-4 -5e-4 -1e-3 -1.5e-3],'LineColor','black');
v = [0.1175 0.115 0.11 0.1 9e-2 7e-2 5e-2 3e-2 1e-2 -1e-4 -5e-4 -1.5e-3];
clabel(C,h,v)
xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')
f = gcf; filename = "Streamlines_Contour_Re_"+Re+"_N_"+N+".png";
file = fullfile("C:\Users\Senne Hemelaar\OneDrive\Documenten\GitHub\CFD_2\Main\Images\Contours", filename);
exportgraphics(f,file,'Resolution',600)
end

%%%====================== Line Plots Constant x ========================%%%
color_string = linspecer(6);

%%% VELOCITY %%%
if (settings.velocity_line_plot_Cx)
figure(4)
hold on; grid on; box on;
if i == 1
plot(ref_Cx.y,ref_Cx.u,'ks-','Linewidth',1.1)
end
plot(results.x(2:end-1),results.u_x05,'color',color_string(i,:),'Linewidth',1.1)
xlabel('$y$','interpreter','latex')
ylabel('$u$','interpreter','latex')
legend('Reference','$N=15$','$N=31$','$N=63$','interpreter','latex','Location','SouthWest')
f = gcf; filename = "Velocity_Reference_x05_Re_"+Re+".png";
file = fullfile("C:\Users\Senne Hemelaar\OneDrive\Documenten\GitHub\CFD_2\Main\Images\Line_Plots", filename);
exportgraphics(f,file,'Resolution',600)
end

%%% PRESSURE %%%
if (settings.pressure_line_plot_Cx)
figure(5)
hold on; grid on; box on;
if i == 1
plot(ref_Cx.y,ref_Cx.p,'ks-','Linewidth',1.1)
end
plot(results.x(2:end-1),results.p_x05,'color',color_string(i,:),'Linewidth',1.1)
xlabel('$y$','interpreter','latex')
ylabel('$p$','interpreter','latex')
ylim([-0.02 0.12])
legend('Reference','$N=15$','$N=31$','$N=63$','interpreter','latex','Location','NorthEast')
f = gcf; filename = "Pressure_Reference_x05_Re_"+Re+".png";
file = fullfile("C:\Users\Senne Hemelaar\OneDrive\Documenten\GitHub\CFD_2\Main\Images\Line_Plots", filename);
exportgraphics(f,file,'Resolution',600)
end

%%% VORTICITY %%%
if (settings.vorticity_line_plot_Cx)
figure(6)
hold on; grid on; box on;
if i == 1
plot(ref_Cx.y,ref_Cx.xi,'ks-','Linewidth',1.1)
end
plot(results.tx,results.xi_x05,'color',color_string(i,:),'Linewidth',1.1)
xlabel('$y$','interpreter','latex')
ylabel('$\xi$','interpreter','latex')
legend('Reference','$N=15$','$N=31$','$N=63$','interpreter','latex','Location','NorthWest')
f = gcf; filename = "Vorticity_Reference_x05_Re_"+Re+".png";
file = fullfile("C:\Users\Senne Hemelaar\OneDrive\Documenten\GitHub\CFD_2\Main\Images\Line_Plots", filename);
exportgraphics(f,file,'Resolution',600)
end

%%%====================== Line Plots Constant y ========================%%%

%%% VELOCITY %%%
if (settings.velocity_line_plot_Cy)
figure(7)
hold on; grid on; box on;
if i == 1
plot(ref_Cy.x,ref_Cy.v,'ks-','Linewidth',1.1)
end
plot(results.x(2:end-1),results.v_y05,'color',color_string(i,:),'Linewidth',1.1)
xlabel('$x$','interpreter','latex')
ylabel('$v$','interpreter','latex')
legend('Reference','$N=15$','$N=31$','$N=63$','interpreter','latex','Location','NorthWest')
f = gcf; filename = "Velocity_Reference_y05_Re_"+Re+".png";
file = fullfile("C:\Users\Senne Hemelaar\OneDrive\Documenten\GitHub\CFD_2\Main\Images\Line_Plots", filename);
exportgraphics(f,file,'Resolution',600)
end

%%% PRESSURE %%%
if (settings.pressure_line_plot_Cy)
figure(8)
hold on; grid on; box on;
if i == 1
plot(ref_Cy.x,ref_Cy.p,'ks-','Linewidth',1.1)
end
plot(results.x(2:end-1),results.p_y05,'color',color_string(i,:),'Linewidth',1.1)
xlabel('$x$','interpreter','latex')
ylabel('$p$','interpreter','latex')
ylim([-0.02 0.12])
legend('Reference','$N=15$','$N=31$','$N=63$','interpreter','latex','Location','North')
f = gcf; filename = "Pressure_Reference_y05_Re_"+Re+".png";
file = fullfile("C:\Users\Senne Hemelaar\OneDrive\Documenten\GitHub\CFD_2\Main\Images\Line_Plots", filename);
exportgraphics(f,file,'Resolution',600)
end

%%% VORTICITY %%%
if (settings.vorticity_line_plot_Cy)
figure(9)
hold on; grid on; box on;
if i == 1
plot(ref_Cy.x,ref_Cy.xi,'ks-','Linewidth',1.1)
end
plot(results.tx,results.xi_y05,'color',color_string(i,:),'Linewidth',1.1)
xlabel('$x$','interpreter','latex')
ylabel('$\xi$','interpreter','latex')
legend('Reference','$N=15$','$N=31$','$N=63$','interpreter','latex','Location','South')
f = gcf; filename = "Vorticity_Reference_y05_Re_"+Re+".png";
file = fullfile("C:\Users\Senne Hemelaar\OneDrive\Documenten\GitHub\CFD_2\Main\Images\Line_Plots", filename);
exportgraphics(f,file,'Resolution',600)
end

end
close all



















