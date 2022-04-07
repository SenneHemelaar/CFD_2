function [E21] = Create_E21(N)
%E21:  Assembles the incidence matrix relating circulation surfaces and line
%      segments

% Horizontal Fluxes
lower_diag = ones((N+1)^2 + N+1,1);
upper_diag = -ones((N+1)^2,1);
u_part = diag(lower_diag) + diag(upper_diag,N+1);
u_part = u_part(1:(N+1)^2,:);

% Vertical fluxes
lower_diag = -ones(N+2,1);
upper_diag = ones(N+1,1);
diag_u     = diag(lower_diag) + diag(upper_diag,1);
diag_u     = diag_u(1:end-1,:);
I          = eye(N+1);
v_part     = kron(I,diag_u);

E21 = [u_part, v_part];

end

