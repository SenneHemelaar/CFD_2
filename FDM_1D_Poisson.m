clear;
close all;
clc;

%   Sample problem, 1D Poisson equaltion
%       d2u/dx2 = f(x)
%       u(0) = 0;   u'(1)=0
%
%   20 NOvember 2017, Marc Gerritsma
%

N = 100;                 % Number of unknowns
x = linspace(0,1,N+1);  % Grid
h = 1/N;                % mesh size

%
%   Setting up the (sparse) matrix A
%

ld = 1/(h*h)*ones(1,N-1);  ld = [ld 0]; ld(N-1) = 2*ld(N-1);
cd = -2/(h*h)*ones(1,N);
td = 1/(h*h)*ones(1,N-1);  td = [0 td];

B=[ld; cd; td];   A = spdiags(B',[-1,0,1],N,N);

%
%   The right hand side vector
%

%   rhs = x(2:N+1)';
rhs = sin(1.5*pi*x(2:N+1)');

sol = A\rhs;

sol = [0; sol];

Nf=100;
xf = linspace(0,1,Nf+1);
%   ex = (xf.^3 - 3*xf)/6;
ex = -4*sin(1.5*pi*xf)/(9*pi^2);

err = -4*sin(1.5*pi*x)/(9*pi^2) - sol';

figure(1)
plot(x,sol,'r*')
hold on
plot(xf,ex)
title('Exact and approximate solution')

figure(2)
plot(x,err,'r*')
title('Error in the points')

L2 = sqrt(err*err'/N)

kin = 0.0;
for i=2:N+1
    kin = kin + 0.5*N*(sol(i)-sol(i-1))^2;
end
kin
   





