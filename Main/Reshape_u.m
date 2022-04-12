function [results, u,v] = Reshape_u(results, N)
%RESHAPE_U Reshapes u vector and turns fluxes into velocity

% Seperate u and v from u vector and reshape
results.v     = results.u(length(results.u)/2+1:end);
results.u     = results.u(1:length(results.u)/2);
results.u     = transpose(reshape(results.u,[N+1,N]));
results.v     = transpose(reshape(results.v,[N,N+1]));

% Devide by respective edges to obtain velocity
for i = 1:length(results.th)
    results.u(i,:) = results.u(i,:)./results.th(i);
    results.v(:,i) = results.v(:,i)./results.th(i);
end

% Take averages to obtain velocities at dual grid nodes
u = zeros(N); v = zeros(N);
for i = 1:N
    u(:,i) = (results.u(:,i)+results.u(:,i+1))./2;
    v(i,:) = (results.v(i,:)+results.v(i+1,:))./2;
end
results.u = u;
results.v = v;

% Create additional vector containing absolute velocity
results.u_abs = zeros(N,N);
for i = 1:N
    for j = 1:N
        results.u_abs(i,j) = sqrt(results.u(i,j)^2 + results.v(i,j)^2);
    end
end

end

