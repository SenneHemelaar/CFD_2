function [results] = Reshape_p(results, N)
%RESHAPE_P Reshapes the pressure vector

results.p = results.p(N+1:end-N);
results.p = reshape(results.p,[N+2,N]);
results.p = transpose(results.p);
results.p = results.p(:,2:end-1);

end

