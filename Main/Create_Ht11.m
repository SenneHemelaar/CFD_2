function [Ht11] = Create_Ht11(N, th, h)
%CREATEHT11 Creates Hodge Matrix Ht11

% Create Diagonal
count = 0;
for i = 1:N
    for j = 1:N+1
        count = count+1;
        Ht11(count) = th(i)/h(j);
    end
end
for j = 1:N+1
    for i = 1:N
        count = count+1;
        Ht11(count) = th(i)/h(j);
    end
end
% Assemble matrix
Ht11 = diag(Ht11);
Ht11 = sparse(Ht11);
end

