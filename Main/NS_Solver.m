function [results, tx, th, x, h] = NS_Solver(N, Re, tol, dt_factor)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                         CFD 2 NS Solver                             %%%
%%%                      Author: Senne Hemelaar                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
% The system that you need to solve will be singular. Matlab gives you a
% warning at each time step. To switch of this warning, remove the comment
% in the next line
%}
warning off
Delta = 1/N;      % uniform spacing to be used in the mapping to compute tx
% wall velocities
BC.U_wall_top   = -1;
BC.U_wall_bot   = 0;
BC.U_wall_left  = 0;
BC.U_wall_right = 0;
BC.V_wall_top   = 0;
BC.V_wall_bot   = 0;
BC.V_wall_left  = 0;
BC.V_wall_right = 0;

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

dt = dt_factor * min(min(h),0.5*Re*min(h)^2);

u = zeros(2*N*(N+1),1);

%%%================ Set up the incindence matrix 'tE21' ================%%%
tE21 = Create_tE21(N);

%%%=================== Seperate 'tE21' and 'u_norm' ====================%%%
u_vec_norm = zeros(4*N,1);
for i = 1:N
    u_vec_norm(1+2*(i-1)) = BC.U_wall_left;
    u_vec_norm(2*i)       = BC.U_wall_right;
    u_vec_norm(2*N+i)     = BC.V_wall_bot;
    u_vec_norm(3*N+i)     = BC.V_wall_top;
end
[tE21, u_norm] = Seperate_tE21(tE21, N, u_vec_norm);

%%%=============== Set up the incidence matrix 'tE10' ==================%%%
[tE10] = Create_tE10(N);

%%%================ Set up the incidence matrix 'E10' ==================%%%
E10 = -tE21';

%%%================ Set up the incidence matrix 'E21' ==================%%%
E21 = Create_E21(N);

%%%=================== Seperate 'E21' and 'u_pres' =====================%%%
u_vec_tan = zeros((N+1)*4,1);
for i = 1:N+1
    u_vec_tan(i)                 = BC.U_wall_bot * h(i);
    u_vec_tan((N+1)+i)           = BC.U_wall_top * h(i);
    u_vec_tan((N+1)*2+1+(i-1)*2) = BC.V_wall_left * h(i);
    u_vec_tan((N+1)*2+(2*i))     = BC.V_wall_right * h(i);
end
[E21, u_pres] = Seperate_E21(E21, N, u_vec_tan);

%%%=============== Set up the Hodge matrices Ht11 and H1t1 =============%%%
Ht11 = Create_Ht11(N, th, h);

H1t1 = zeros(2*N*(N+1), 2*N*(N+1));
for i = 1:length(Ht11)
    H1t1(i,i) = 1/Ht11(i,i);
end
H1t1 = sparse(H1t1);

%%%===================== Set up the Hodge matrix Ht02 ==================%%%
[Ht02] = Create_Ht02(N, h);

H2t0 = zeros(length(Ht02),length(Ht02));
for i = 1:length(Ht02)
    H2t0(i,i) = 1/Ht02(i,i);
end
H1t1 = sparse(H1t1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%=================== Non adjustable part of the code =================%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u_pres_vort=Ht02*u_pres; %U_pres to outer oriented 0 form representing contribution of boundary conditions to point wise vorticity
u_pres = H1t1*E21'*Ht02*u_pres; %U_pres to inner oriented 1 forms

%{
% Now all matrices are set up and the time stepping can start. 'iter' will
% record the number of time steps. This allows you to give output after a
% preselected number of time steps.
% 
% 'diff' will be the maximal du/dt or dv/dt. If 'diff' is sufficiently
% small, steady state has been reached. Determine a suitable value for
% 'tol'
%} 

convective = zeros(2*N*(N+1),1);
ux_xi = zeros((N+1)*(N+1),1);
uy_xi = zeros((N+1)*(N+1),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%======================== Time stepping loop =========================%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
diff = 1;
iter = 1;
m    = 1;

%%%============ Set up the matrix for the Poisson equation =============%%%  
A = -tE21*Ht11*tE21';

%%%====== Perform an LU-decomposition for the pressure matrix A=========%%%
[L,U] = lu(A);

%%%== Abbrev. for some matrix products that are constant in the loop ===%%%
VLaplace = H1t1*E21'*Ht02*E21;
DIV = tE21*Ht11;

while diff > tol
%while iter < 10
     %{   
    %Vector xi is obtained. It corresponds with the point-wise vorticity
    %at each cell
    
    %Vectors ux_xi and uy_xi correspond with the multiplications of
    %xi with the horizontal and vertical velocity components at each cell.
    %Only the cells required for vector convective are calculated. The
    %ordering of each vector with respect to the ordering of cells in the
    %grid is different (left to right for ux_xi and bottom to top for
    %uy_xi)
    %}

    xi = Ht02*E21*u + u_pres_vort;
    
    for i=1:N+1
        for j=1:N+1
            k = j + (i-1)*(N+1); 
            if j==1
                ux_xi(k) = BC.U_wall_bot*xi(i+(j-1)*(N+1));    
                uy_xi(k) = BC.V_wall_left*xi(j+(i-1)*(N+1));
            elseif j==N+1
                ux_xi(k) = BC.U_wall_top*xi(i+(j-1)*(N+1));
                uy_xi(k) = BC.V_wall_right*xi(j+(i-1)*(N+1));
            else
                ux_xi(k) = (u(i+(j-1)*(N+1))+u(i+(j-2)*(N+1)))*xi(i+(j-1)*(N+1))/(2.*h(i));                      
                uy_xi(k) = (u(N*(N+1)+j+(i-1)*N) + u(N*(N+1)+j-1+(i-1)*N))*xi(j+(i-1)*(N+1))/(2.*h(i));  
            end
        end
    end

    for  i=1:N
        for j=1:N+1
            convective(j+(i-1)*(N+1)) = -(uy_xi(j+(i-1)*(N+1))+uy_xi(j+i*(N+1)))*h(j)/2.;
            convective(N*(N+1)+i+(j-1)*N) = (ux_xi(j+(i-1)*(N+1))+ux_xi(j+i*(N+1)))*h(j)/2.;
        end
    end
    

%%%====== Set up the RHS for the Poisson equation for the pressure =====%%%
    rhs_Poisson = DIV*(u/dt  - convective - VLaplace*u/Re - u_pres/Re)...
                  + u_norm/dt;  

%%%===================== Solve for the new pressure ====================%%%
    temp = L\rhs_Poisson;
    p = U\temp;
    
%%%==== Store the velocity from prev. time step in the vector u_old ====%%%
    uold = u;
    
%%%====================== Udate the velocity field =====================%%%    
    u = u - dt* (convective - tE21'*p + VLaplace*u/Re + u_pres/Re); 
    
    %{
    %  Every other 1000 iterations check whether you approach steady state
    %  and check whether yopu satisfy conservation of mass. The largest
    %  rate at which mass is destroyed or created is denoted by 'maxdiv'.
    %  This number should be very small, in the order of machine precision.
    %}
    if mod(iter,1000) == 0
    
        maxdiv = max(DIV*u + u_norm); 
        
        diff = max(abs(u-uold))/dt
        
        results.t(m) = toc;
        
        m = m + 1;
    end
    diff_list(iter) = diff;
    iter = iter + 1;
end

results.u     = u;
results.p     = p;
results.xi    = xi;
results.ux_xi = ux_xi;
results.uy_xi = uy_xi;
results.N     = N;
results.Re    = Re;
results.tx    = tx;
results.th    = th;
results.x     = x;
results.h     = h;
results.diff  = diff_list;
results.iter  = iter;
results.dt    = dt;
results.tE10  = tE10;
results.tE21  = tE21;
results.E10   = E10;
results.E21   = E21;
results.Ht11  = Ht11;
results.Ht02  = Ht02;
results.H2t0  = H2t0;
results.A     = A;
end



