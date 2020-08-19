% script name: "test_first_moment_versions"
%

% load coeffs and gamma 
if exist('small_example.mat')
    load('small_example.mat');
else
    prepare_coef_and_gamma_example;
    clear;
    load('small_example.mat');
end

% random distribution
P = 4;
B = cell(P,1);
B{1} = 1;

for j=2:P
    B{j} = randn(2*j-1,2*j-1);
end

K = size(gamma.coeff{1},2); %size(A{1}{1},1);

% compare
tic
[mu1] = FirstMoment_PSWF_naive_fixed2(A, B, gamma);
toc();

tic
L = max(gamma.band_idx_3d);
[sign_mat, Gamma_mat] = preprocessing_mu1_coef(P, L, gamma);
% [mu1_new] = FirstMoment_PSWF_v1(A, B, gamma);% Gamma_mat, sign_mat);
[mu1_new] = FirstMoment_PSWF_v2(A, B, Gamma_mat, sign_mat);

toc();


norm(mu1(:)-mu1_new(:))

