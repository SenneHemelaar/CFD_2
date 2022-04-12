%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                           CFD 2 Plotter                             %%%
%%%                      Author: Senne Hemelaar                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; 
% close all; clc

%%%========================== Load Results =============================%%%
N = 55;
Re = 1000;
matfile = "Results_N_"+N+"_Re_"+Re+".mat";
load(matfile)

%%%======================== Plotting Routines ===========================%%%
[x,y] = meshgrid(results.x(2:end-1),results.x(2:end-1));
figure(1)
contourf(x,y,results.p,[-0.002 0.0 0.02 0.05 0.07 0.09 0.11 0.12 0.17 0.3],...
         'FaceColor','white','ShowText','on')
xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')
title('$p$','interpreter','latex')

% [x,y] = meshgrid(results.tx,results.tx);
% figure(2)
% contourf(x,y,results.xi,[-3.0 -2.0 -1.0 -0.5 0.0 0.5 1.0 2.0 3.0 4.0 5.0]...
%          ,'FaceColor','white','ShowText','on')
% xlabel('$x$','interpreter','latex')
% ylabel('$y$','interpreter','latex')
% title('$\omega$','interpreter','latex')

% [x,y] = meshgrid(results.x(2:end-1),results.x(2:end-1));
% figure(4)
% % contourf(x,y,results.u_abs,[0.001 0.01 0.02 0.05 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1],...
% %          'FaceColor','white','ShowText','on')
% streamline(x,y,results.u,results.v,x,y)
% xlabel('$x$','interpreter','latex')
% ylabel('$y$','interpreter','latex')
% title('$u$','interpreter','latex')






