% preliminaries
fprintf('\n  \n  \n  \n')
clear;

%% generating "ground truth"
fprintf('**** Starting simulation ****\n');
fprintf('\n  \n')
fprintf('loading data..')
P  = 3;     % aribitrary, chosen to be small
B = get_B_inplane(P);
% bval1 = draw_inplane_dist(B,1, -50, -70, 'inplane_dist_most_nu');

% B = totally_NonUniform_SOS_dist(P);%     % "Project" to positive
% 
% % getting random inplane uniform distributon -- expansion coefficients
% % bval1 = draw_inplane_dist(B,1, -50, -70, 'inplane_dist_most_nu');
load ('SO3_fifteen.mat');
[AI,AE,Acon] = linear_cons_B2(numel(B)-1,SO3);
[BB,Breal]    = project_B_to_positive2(AI,AE,Acon,B);
scl = 1/B{1};
for i=1:numel(B)
    B{i} =  scl*B{i};
end
%% Simulate kspace function
gridSize = 19;  % number of voxels in each dimenison. take odd.

% Sum of Gaussian
% vol1    = cryo_gaussian_phantom_3d('C1_params',gridSize,1);
sigma = 200;
% T     = (gridSize/5)*[0.1 -0.2 -0.4; 0.0 0.2 0; 0.1 -0.4 -0.45]';%[0.1 -0.2 -0.4; 0.0 0.2 0; 0.1 -0.4 -0.45]';% ; 0.13 -0.2 0.1;0.1 0 -0.15]';
% g     =  @(k) exp(1i*k'*T(:,1)).*exp(-pi^2/sigma*sum(k.*k)).'  + ...
%     exp(1i*k'*T(:,2)).*exp(-pi^2/sigma*sum(k.*k)).'  + ...
%     exp(1i*k'*T(:,3)).*exp(-pi^2/sigma*sum(k.*k)).' ;

T     = (gridSize/5)*[0, 0, 0;
                      0.15, 0.10, 0.10;
                      -0.1, -0.15, 0;
                      -0.1, 0.05, -0.05;
                      0, 0.1, -0.1;]';%[0.1 -0.2 -0.4; 0.0 0.2 0; 0.1 -0.4 -0.45]';% ; 0.13 -0.2 0.1;0.1 0 -0.15]';
g     =  @(k) exp(1 * 1i*k'*T(:,1)) .* exp(-0.7   * pi^2/sigma*sum(k.*k)).' + ...
              exp(1 * 1i*k'*T(:,2)) .* exp(-0.5 * pi^2/sigma*sum(k.*k)).' + ...
              exp(1 * 1i*k'*T(:,3)) .* exp(-0.4   * pi^2/sigma*sum(k.*k)).' + ...
              exp(1 * 1i*k'*T(:,4)) .* exp(-0.6   * pi^2/sigma*sum(k.*k)).' + ...
              exp(1 * 1i*k'*T(:,5)) .* exp(-0.8   * pi^2/sigma*sum(k.*k)).' ;

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

%% calculate 3D expansion and gamma coefficients
beta  = .85;       % Bandlimit ratio
delta = 0.99;    % Truncation parameter
eps_p = 1e-3;    % Prescribed accuracy

fprintf('Calculating 3D coefficients...');
AFull = pswf_t_f_3d(vol, beta, delta);
fprintf('DONE \n');

fprintf('Calculating gamma coefficients...');
radius = floor(gridSize/2);
c     = beta*pi*radius;              % Nyquist bandlimit
gamma = PSWF_2D_3D_T_mat(c, delta, eps_p);
% Trncated
[gamma,A] = gamma_truncate_2(gamma,A);

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

fprintf('Moments calculation...');
% first moment
[ m1_true ] = FirstMoment_PSWF_v2(A, B, Gamma_mat, sign_mat);
fprintf('First moment is DONE. \n');

% second moment
fprintf('Mu2 calculation...');
tic
[m2_true] = SecondMoment_PSWF_v2(A, B, gamma, C_tensor, Gamma_mat);
tt = toc;
fprintf('DONE in about %d seconds \n', round(tt));


%% Compute Projecitons
paramsName = 'C1_params' ;
Ns = floor(logspace(3,6));
repNumber = 10;
errorM1 = zeros(length(Ns), repNumber);
errorM2 = zeros(length(Ns), repNumber);

for indexN = 1 : length(Ns)
    total_N = Ns(indexN);
    tmpErrorM1 = zeros(length(Ns),1);
    tmpErrorM2 = zeros(length(Ns),1);
    for iRep = 1 : repNumber

        [~, proj_PSWF, weight] = GenerateObservationsPSWF(paramsName,...
                                        total_N, gridSize, B, beta, eps_p, gamma, g, x_2d, y_2d);
        %% Estimate m1_hat and m2_hat
        [m1_hat, m2_hat] = EstimateMoments(proj_PSWF, weight, total_N, gamma);
        tmpErrorM1(iRep) = norm(m1_hat(:) - m1_true(:)).^2;
        tmpErrorM2(iRep) = norm(m2_hat(:) - m2_true(:)).^2;
    end
    errorM1(indexN,:) = tmpErrorM1;
    errorM2(indexN,:) = tmpErrorM2;

end

%% Display
figure;
subplot(2,1,1);
loglog(Ns, mean(errorM1,2),'*-');
title('Convargance for M1');

subplot(2,1,2);
loglog(Ns, mean(errorM2,2),'*-');
title('Convargance for M2');
