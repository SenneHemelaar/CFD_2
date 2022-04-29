function [results] = Reshape_p(results, N)
%RESHAPE_P Reshapes the pressure vector

results.p = results.p(N+1:end-N);
results.p = reshape(results.p,[N+2,N]);
results.p = transpose(results.p);
results.p = results.p(:,2:end-1);

% Convert from total to static pressure
for i =1:N
    for j = 1:N
        results.p(i,j) = results.p(i,j) - 0.5*results.abs_u(i,j)^2;
    end
end

% Find the center value for p and substract from p matrix
p_center = results.p(round(N/2),round(N/2));
results.p = results.p - p_center;
end

