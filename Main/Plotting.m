%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                           CFD 2 Plotter                             %%%
%%%                      Author: Senne Hemelaar                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc

%%%========================== Load Results =============================%%%
N = 50;
Delta = 1/N;  

u  = load('u_N_50');  u  = u.u;
p  = load('p_N_50');  p  = p.p;
vort = load('xi_N_50'); vort = vort.xi;
tx = zeros(1,N+1);
for i=1:N+1
    xi = (i-1)*Delta;
    tx(i) = 0.5*(1. - cos(pi*xi));
end
th = zeros(N,1);
th = tx(2:N+1) - tx(1:N);
x1 = 0.5*(tx(1:N) + tx(2:N+1));

% Vorticity
vort = transpose(reshape(vort,[51,51]));

% pressure
p_bot = [0, p(1:50)', 0];
p_top = [0, p(end-49:end)', 0];
p_mid = p(51:end-50);
p_mid = transpose(reshape(p_mid,[52,50]));
p = [p_top; p_mid; p_bot];
x2 = 0.5*(tx(1:N) + tx(2:N+1));
x2 = [0 x2 1];
[x,y] = meshgrid(x2,x2);
figure(3)
contourf(x,y,p,[-0.002 0.0 0.02 0.05 0.07 0.09 0.11 0.12 0.17 0.3],...
         'FaceColor','white','ShowText','on')
xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')
title('$p$','interpreter','latex')

% velocity
v = u(length(u)/2+1:end);
u = u(1:length(u)/2);
u = transpose(reshape(u,[51,50]));
v = transpose(reshape(v,[50,51]));

% [x,y] = meshgrid(tx,x1);
% figure(1)
% contourf(x,y,u)
% colorbar
% xlabel('$x$','interpreter','latex')
% ylabel('$y$','interpreter','latex')
% title('$u$','interpreter','latex')
% 
% [x,y] = meshgrid(x1,tx);
% figure(2)
% contourf(x,y,v)
% colorbar
% xlabel('$x$','interpreter','latex')
% ylabel('$y$','interpreter','latex')
% title('$v$','interpreter','latex')


% [x,y] = meshgrid(tx, tx);
% figure(4)
% contourf(x,y,vort,[-3.0 -2.0 -1.0 -0.5 0.0 0.5 1.0 2.0 3.0 4.0 5.0],'FaceColor','white','ShowText','on')
% xlabel('$x$','interpreter','latex')
% ylabel('$y$','interpreter','latex')
% title('$\omega$','interpreter','latex')





