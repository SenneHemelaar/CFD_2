function [Ht02] = Create_Ht02(N, h)
%CREATEHT11 Creates Hodge Matrix Ht02

Ht02 = zeros((N+1)^2, (N+1)^2);
for i = 1:(N+1)
    for j = 1:(N+1)
        Ht02((N+1)*(i-1)+(j-1) + 1, (N+1)*(i-1)+(j-1)+1)...
        = 1 / (h(i) * h(j));
    end
end
Ht02 = sparse(Ht02);
end

