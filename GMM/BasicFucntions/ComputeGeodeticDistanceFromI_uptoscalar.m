function [delta2] = ComputeGeodeticDistanceFromI_uptoscalar(A)
% A and B same szie square, semi-postaive definte, stmetric matrices
A = A / norm(A, 'fro') * sqrt(size(A,1));
egis = eig(inv(A));
% eigsInv = 1 ./ egis;
delta2 = sqrt(sum(log(egis).^2));
% geoDistTmp = sqrt(sum(log(egis).^2));
% % geoDistTmp = sqrt(size(A,1)) * log(2);
% B = A * exp(-1 * geoDistTmp  / norm(A, 'fro'));% (sqrt(size(A,1))));% * norm(A, 'fro')));
% % B = B / norm(B, 'fro');
% eigB = eig(inv(B));
% % eigBInv = 1 ./ eigB;
% delta2 = sqrt(sum(log(eigB).^2));

end

