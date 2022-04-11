%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                           CFD 2 Plotter                             %%%
%%%                      Author: Senne Hemelaar                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc

%%%========================== Load Results =============================%%%
load('Results_N_63_Re_1000')

% [x,y] = meshgrid(results.x(2:end-1),results.x(2:end-1));
% figure(1)
% contourf(x,y,results.p,[-0.14 -0.12 -0.1 -0.08 -0.06 -0.03 -0.027 -0.024 -0.2 -0.017 -0.014 -0.01 -0.007 -0.004 0.0 0.02 0.05 0.07 0.09 0.11 0.12 0.17 0.3],...
%          'FaceColor','white','ShowText','on')
% xlabel('$x$','interpreter','latex')
% ylabel('$y$','interpreter','latex')
% title('$p$','interpreter','latex')

[x,y] = meshgrid(results.tx,results.tx);
figure(2)
contourf(x,y,results.xi,[-3.0 -2.0 -1.0 -0.5 0.0 0.5 1.0 2.0 3.0 4.0 5.0]...
         ,'FaceColor','white','ShowText','on')
xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')
title('$\omega$','interpreter','latex')

