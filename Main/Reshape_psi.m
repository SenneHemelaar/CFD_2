function [results] = Reshape_psi(results, N)
%RESHAPE_P Reshapes the stream function

psi = linsolve(results.tE10,results.tu);
results.psi = transpose(reshape(psi,[N+1,N+1])); 

end