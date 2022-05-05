clear all; close all;


N_list = 20:1:30;
Re = 1000;

for j = 1:length(N_list)
    N = N_list(j);
    Delta = 1/N;
tx = zeros(1,N+1);
for i=1:N+1
    xi = (i-1)*Delta;
    tx(i) = 0.5*(1. - cos(pi*xi));
end
th = zeros(N,1);
th = tx(2:N+1) - tx(1:N);
x = 0.5*(tx(1:N) + tx(2:N+1));
x = [0 x 1];
h = zeros(N+1,1);
h = x(2:N+2) - x(1:N+1);

min_h(j) = min(h);
end
Re_part = 0.5*Re*min_h.^2;

figure;
hold on; grid on; box on;
plot(N_list,min_h,'k','linewidth',1.1)
plot(N_list,0.5*Re*min_h.^2,'k--','linewidth',1.1)
xlabel('N','interpreter','latex')
ylabel('$\Delta t$ [s]','interpreter','latex')
legend('$h_{min}$','$0.5Re\,h_{min}^{2}$','interpreter','latex','Fontsize',11)
set(gcf,'position',[100,100,750,400])





