clear; close all

N = 4;
S = N*N + 4*N;

% Horizontal fluxes inner grid
upper_diag = ones(N,1);
lower_diag = -ones(N+1,1);
diag_u = diag(lower_diag) + diag(upper_diag,1);
diag_u = diag_u(1:end-1,:);
I = eye(N);
u_part_inner = kron(I,diag_u);

% Horizontal fluxes virtual surfaces
u_part_outer = zeros(N*N+N, N*N+N);
u_part_outer(1,1) = 1;
u_part_outer(2*N,end) = -1;
for i = 1:N-1
   u_part_outer(2+(N-1)*(i-1),(N+1)*i) = -1;
   u_part_outer(3+(N-1)*(i-1),(N+1)*i+1) = 1;
end

% Vertical fluxes inner grid
lower_diag = -ones(N*N + N,1);
upper_diag = ones(N*N,1);
v_part_inner = diag(lower_diag) + diag(upper_diag,N);
v_part_inner = v_part_inner(1:N*N,:);

% Vertical fluxes virtual surfaces
v_part_outer = zeros(2*N,N*N+N);
for i = 1:N
    v_part_outer(i,i) = 1;
    v_part_outer(end-(i-1),end-(i-1)) = -1;
end
v_part_outer = [zeros(size(v_part_outer)) ; v_part_outer];

% Assemble tE21 matrix
% u = [u_part_inner ; u_part_outer];
% v = [v_part_inner ; v_part_outer];
% tE21 = [u, v];







% Vertical fluxes
% lower_diag = -ones(N*(N+3),1);
% upper_diag = ones(N*(N+3) - N,1);
% diag_v     = diag(lower_diag) + diag(upper_diag,N);
% diag_v     = diag_v(1:end-N,:);
% v          = diag_v;
% v = [v(1:N,:); zeros(1,length(v)); v(N+1:end,:)];
% v = [v(1:end-N,:); zeros(1,length(v)); v(end-N+1:end,:)];
% for i = 1:N-1
% v = [v(1:2*N+1 + (i-1)*(N+2),:); zeros(2,N*(N+3)); v(2*N+2+(i-1)*(N+2):end,:)];
% end
% 
% % Assemble tE21 matrix
% tE21 = [u, v];
% % tE21 = sparse(tE21);
% 
% tE21_noboundary = [tE21(:,2:5) , tE21(:,8:11) , tE21(:,14:17) , tE21(:,22:33) ];
% 
% % E21 = [tE21(1:N+1,:) ; tE21(2*N+2:2*N+3,:) ; tE21(3*N+4:3*N+5,:) ; tE21(end-N:end,:)];
% % E21 = [tE21(N+2:2*N+1,:) ; tE21(2*N+4:3*N+3,:) ; tE21(3*N+6:end-N-1,:)];





