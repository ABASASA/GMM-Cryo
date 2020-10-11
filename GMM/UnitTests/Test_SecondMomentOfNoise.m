% preliminaries
fprintf('\n  \n  \n  \n')
clear;

%% generating "ground truth"
fprintf('**** Starting simulation ****\n');
fprintf('\n  \n')
fprintf('loading data..')
P  = 3;     % aribitrary, chosen to be small
SNR = 5;
B = get_B_inplane(P);
Ns = 1;
beta  = 1;       % Bandlimit ratio
delta = 0.99;    % Truncation parameter
eps_p = 1e-3;    % Prescribed accuracy

%% Simulate kspace function
gridSize = 15;  % number of voxels in each dimenison. take odd.

sigma = 200;

T     = (gridSize/5)*[0 0 0; 0.08 0.1 0; -0.1 0 0.1]';% ; 0.13 -0.2 0.1;0.1 0 -0.15]';
g     =  @(k) exp(1i*k'*T(:,1)).*exp(-pi^2/sigma*sum(k.*k)).'  + ...
    exp(1i*k'*T(:,2)).*exp(-pi^2/sigma*sum(k.*k)).'  + ...
    exp(1i*k'*T(:,3)).*exp(-pi^2/sigma*sum(k.*k)).' ;%+ ...
%% Vol on kspace over cartesian grid

radius   = floor(gridSize/2);
if mod(gridSize,2)==0
    x_1d = (-radius:1:radius-1);%/radius;   % - Even number of points
else
    x_1d = (-radius:1:radius);%/radius;   % - Odd number of points
end

% the 3d grid
[x_3d,y_3d,z_3d] = meshgrid(x_1d,x_1d,x_1d);
vec_3d_grid      = [x_3d(:),y_3d(:),z_3d(:)].'*pi/2;

% parameterization of the 2d slices
x_2d = x_3d(:,:,1);
y_2d = y_3d(:,:,1);

% evaluate the function
volf = g(vec_3d_grid);

% back to real domain
volk = reshape(volf, gridSize, gridSize, gridSize);
vol  = real(fftshift(ifftn(ifftshift(volk)))); % in real domain
vol  = vol*floor(size(x_2d,1)/2);   % Nir Sunday ?

%% calculate 3D expansion and gamma coefficients


fprintf('Calculating 3D coefficients...');
AFull = pswf_t_f_3d(vol, beta, delta);
fprintf('DONE \n');

fprintf('Calculating gamma coefficients...');
radius = floor(gridSize/2);
c     = beta*pi*radius;              % Nyquist bandlimit
gamma = PSWF_2D_3D_T_mat(c, delta, eps_p);
%Trncated
% [gamma,A] = gamma_truncate_2(gamma,AFull);
A = AFull;
fprintf('DONE \n');

%% calculating moments
fprintf('Calculating moments: \n');

% preprocessing
L = max(gamma.band_idx_3d);  % the overall degree. Common to both 2D and 3D
P = size(B,1);               % expansion length of the distribution
M = max(gamma.ang_idx_2d)+1; % angular size of the moment

fprintf('Run preprocessing...');
tic
[sign_mat, Gamma_mat] = preprocessing_mu1_coef(P, L, gamma);
[C_tensor]            = preprocessing_mu2_coef_V2(gamma, P, M);
tt = toc;
fprintf('DONE in about %d seconds \n', round(tt));

%% Compute Projecitons
paramsName = [];%'C1_params' ;


total_N = 200000;
[projs ,weight] = ComputeProjection (total_N, gridSize, paramsName, B, g, x_2d, y_2d);
[projsNoised, noise, ~, sigma] = cryo_addnoise(projs, SNR,'gaussian');
sigma = 1;
noise = randn(size(noise)) * sigma;
[proj_PSWFNoiseOnly] = Projection2PSWF(noise, beta, eps_p, gamma);
%%
m2_hat = 0;
for i=1:total_N
    % averaging
    current_vec_proj = proj_PSWFNoiseOnly(:,:,i);    
    m2_hat = m2_hat + current_vec_proj(:) * current_vec_proj(:)' * weight(i);
end

m2Emp = zeros(max(gamma.ang_idx_2d)+1, nnz(gamma.ang_idx_2d==0),...
                            max(gamma.ang_idx_2d)+1, nnz(gamma.ang_idx_2d==0));
for ii = 1 : total_N
    for m=0:max(gamma.ang_idx_2d)
        for k=0:nnz(gamma.ang_idx_2d==0)-1
            m2Emp(m+1,k+1,:,:) = reshape(m2_hat(m+1 + (k)*(max(gamma.ang_idx_2d)+1),:),...
                 max(gamma.ang_idx_2d)+1, nnz(gamma.ang_idx_2d==0) );
        end
    end
end
m2Emp = m2Emp / total_N;
m2Emp1 = m2Emp / sum(weight);