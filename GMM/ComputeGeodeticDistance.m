function [delta2] = ComputeGeodeticDistance(A,B)
% A and B same szie square, semi-postaive definte, stmetric matrices

C = inv(A) * B;

eigsLog2 = log(eig(C)).^2;

delta2 = sum(eigsLog2)^0.5;
end

