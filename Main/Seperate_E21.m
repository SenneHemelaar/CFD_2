function [E21, u_pres] = Seperate_E21(E21, N, u_vec_tan)
%SEPERATE_E21 Seperates the inner fluxes of tE21 and the part where the
% boundary conditions are imposed on.

% Split E21 in horizontal and vertical flux parts
u_part = E21(:,1:(N+1)^2+(N+1));
v_part = E21(:,(N+1)^2+(N+1)+1:end);

% u_pres part
u_pres_u_part = [u_part(:,1:N+1), u_part(:,end-(N):end)];
u_pres_v_part = v_part(:,1);
for i = 1:N
    u_pres_v_part = [u_pres_v_part, v_part(:,(N+2)*i:(N+2)*i+1)];
end
u_pres_v_part = [u_pres_v_part, v_part(:,end)];
u_pres = [u_pres_u_part, u_pres_v_part] * u_vec_tan;

% E21 part
E21_u_part = u_part(:,N+2:end-(N+1));
E21_v_part = v_part(:,2:N+1);
for i = 1:N
    E21_v_part = [E21_v_part, v_part(:,(N+4)+(i-1)*(N+2):(N+4)+(i-1)*(N+2)+(N-1))];
end
E21 = [E21_u_part, E21_v_part];
E21 = sparse(E21);

end

