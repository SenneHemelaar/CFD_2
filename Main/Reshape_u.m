function [results, u,v] = Reshape_u(results, N)
%RESHAPE_U Reshapes u vector and turns fluxes into velocity

% Seperate u and v from u vector and reshape
results.vec_v     = results.tu(length(results.tu)/2+1:end);
results.vec_u     = results.tu(1:length(results.tu)/2);

results.vec_u     = transpose(reshape(results.vec_u,[N+1,N]));
results.vec_v     = transpose(reshape(results.vec_v,[N,N+1]));

% Devide by respective edges to obtain velocity
for i = 1:length(results.th)
    results.vec_u(i,:) = results.vec_u(i,:)./results.th(i);
    results.vec_v(:,i) = results.vec_v(:,i)./results.th(i);
end

% Take averages to obtain velocities at dual grid nodes
u = zeros(N); v = zeros(N);
for i = 1:N
    u(:,i) = (results.vec_u(:,i)+results.vec_u(:,i+1))./2;
    v(i,:) = (results.vec_v(i,:)+results.vec_v(i+1,:))./2;
end
results.vec_u = u;
results.vec_v = v;

% Create additional vector containing absolute velocity
results.u_abs = zeros(N,N);
for i = 1:N
    for j = 1:N
        results.abs_u(i,j) = sqrt(results.vec_u(i,j)^2 + results.vec_v(i,j)^2);
    end
end

end

