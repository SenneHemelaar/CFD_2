clear; close all;

N = 4;

% horizontal part
diag_center = -ones((N+1)*(N+1),1);
diag_upper  = ones(N*(N+1),1);
u_part = diag(diag_center) + diag(diag_upper,N+1);
u_part = u_part(1:N*(N+1),:);

% vertical part
diag_lower = ones(N+1,1);
diag_upper = -ones(N,1);
v_part = diag(diag_lower) + diag(diag_upper,1);
v_part = v_part(1:N,:);
I = eye(N+1);
v_part = kron(I,v_part);

% Assemble matrix
tE10 = [u_part; v_part];


