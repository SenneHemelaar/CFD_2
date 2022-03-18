function [tE21, u_norm] = Seperate_tE21(tE21, N, u_vec)
%SEPERATE_TE21 Seperates the inner fluxes of tE21 and the part where the
% boundary conditions are subjected on.

norm = tE21;
for i = 1:N
    norm = [norm(:,1:1+2*(i-1)), norm(:,N+3 + (2*(i-1)):end)];
end
norm = [norm(:,1:3*N), norm(:,end-(N-1):end)];

u_norm = norm*u_vec;

tE21 = [tE21(:,2:(N+3)*N-1), tE21(:,(N+3)*N+1:end)];
for i = 1:N-1
   tE21 =  [tE21(:,1:N+1 + (i-1)*(N+1)) , tE21(:,N+4 + (i-1)*(N+1):end)];
end
tE21 = tE21(:,1:end-N);
tE21 = [tE21(:,1:N*(N+1)) , tE21(:,N*(N+1) + N+1:end)];

end

