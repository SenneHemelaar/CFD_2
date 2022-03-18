function [tE21] = Create_tE21(N)
%tE21: Assembles the incidence matrix relating the surfaces and line
%      segments

% Horizontal fluxes
upper_diag = ones(N+2,1);
lower_diag = -ones(N+3,1);
diag_u     = diag(lower_diag) + diag(upper_diag,1);
diag_u     = diag_u(1:end-1,:);
I          = eye(N);
u_part     = kron(I,diag_u);
zero_part  = zeros(N,length(u_part));
u          = [zero_part; u_part; zero_part];

% Vertical fluxes
lower_diag = -ones(N*(N+3),1);
upper_diag = ones(N*(N+3) - N,1);
diag_v     = diag(lower_diag) + diag(upper_diag,N);
diag_v     = diag_v(1:end-N,:);
v          = diag_v;
v          = [v(1:N,:); zeros(1,length(v)); v(N+1:end,:)];
v          = [v(1:end-N,:); zeros(1,length(v)); v(end-N+1:end,:)];
for i      = 1:N-1
v          = [v(1:2*N+1 + (i-1)*(N+2),:); zeros(2,N*(N+3));...
              v(2*N+2+(i-1)*(N+2):end,:)];
end

% Assemble tE21 matrix
tE21 = [u, v];
% tE21 = sparse(tE21);

end

