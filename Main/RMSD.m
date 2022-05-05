function [RMSD] = RMSD(y_hat,y,N)

RMSD = zeros(length(y_hat),1);
for i = 1:length(y_hat)
RMSD(i) = (y_hat(i) - y(i))^2;
end
RMSD = sum(RMSD);
RMSD = sqrt(RMSD);
RMSD = RMSD/N^2;

end

