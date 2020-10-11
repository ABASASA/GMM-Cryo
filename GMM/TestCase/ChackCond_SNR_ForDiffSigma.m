% preliminaries
fprintf('\n  \n  \n  \n')
clear;

%% generating "ground truth"
fprintf('**** Starting simulation ****\n');
fprintf('\n  \n')
fprintf('loading data..')

P = 3;
beta  = 1;       % Bandlimit ratio
delta = 0.99;    % Truncation parameter
eps_p = 1e-3;    % Prescribed accuracy
gridSize = 23;  % number of voxels in each dimenison. take odd.
total_N      = 20000; % Number of projection
isPadding = true;
duplicateExpFlag = true;

sigmaScalars = logspace(-3,1,15);
svPercent = inf;

%% Simulate kspace function

% Sum of Gaussian
sigma = 200;


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
vol  = vol*floor(size(x_2d,1)/2);   % Nir Sunday ?



%% calculate 3D expansion and gamma coefficients


fprintf('Calculating 3D coefficients...');
AFull = pswf_t_f_3d(vol, beta, delta);
volBack      = pswf_t_b_3d(AFull,     gridSize, beta, delta);


fprintf('DONE \n');

fprintf('Calculating gamma coefficients...');
radius = floor(gridSize/2);
c     = beta*pi*radius;              % Nyquist bandlimit
gamma = PSWF_2D_3D_T_mat(c, delta, eps_p);

% Trncated
if isPadding
    [gamma,A] = gamma_truncate_2(gamma,AFull);
else
    A = AFull;
end

%% calculating moments
fprintf('Calculating moments: \n');

% preprocessing
L = max(gamma.band_idx_3d);  % the overall degree. Common to both 2D and 3D
M = max(gamma.ang_idx_2d)+1; % angular size of the moment

fprintf('Run preprocessing...');
tic
[sign_mat, Gamma_mat] = preprocessing_mu1_coef(P, L, gamma);
[C_tensor]            = preprocessing_mu2_coef_V2(gamma, P, M);
tt = toc;
fprintf('DONE in about %d seconds \n', round(tt));


%% naming
nameit = ['balls_example_grid_size_',num2str(gridSize)];
tt     = datetime('now');
nameit = [nameit,'P_',num2str(P),'_',num2str(tt.Hour),'_',num2str(tt.Minute)];

%% distribution

load ('SO3_fifteen.mat');

B = get_B_inplane(P);

[AI,AE,Acon] = linear_cons_B2(numel(B)-1,SO3);
[BB,Breal]    = project_B_to_positive2(AI,AE,Acon,B);
scl = 1/B{1};
for i=1:numel(B)
    B{i} =  scl*B{i};
end

%% Compute Projecitons
paramsName = '' ;
conds = zeros(length(sigmaScalars), 1);
SNRs = zeros(length(sigmaScalars), 1);
geoDist = zeros(length(sigmaScalars), 1);
for iSigma = 1 : length(sigmaScalars)
    sigmaScalar = sigmaScalars(iSigma);
%     [sigmaNoise] = RadiallogNoise(gridSize, sigmaScalar);
    [sigmaNoise] = WhiteNoise(gridSize, sigmaScalar);


    [~, proj_PSWF, weight, ~, signalAmp4SNR] = GenerateObservationsPSWF(paramsName,...
                                    total_N, gridSize, B, beta, eps_p, gamma, g, x_2d, y_2d, sigmaNoise);

    %     [Bias] = ComputeBiasMatrixEmpricly(B, sigmaMoment, gridSize,  g, x_2d, y_2d, gamma, BiasNumberOfProjection, numRepBias);
    [proj ,~] = ComputeProjection (1, gridSize, paramsName, B, g, x_2d, y_2d);

    [Bias] = ComputeBiasSecondMomentAnalyticly(size(proj,1), floor(size(proj,1)/2), beta,...
                                                eps_p, [], sigmaNoise, gamma);

    [m1_hat, m2_hat] = EstimateMoments(proj_PSWF, weight, total_N, gamma, Bias);

    %% Create trimming data
    m1Size = size(m1_hat);
    m2Size = size(m2_hat);

    [vecBoolMoments, ~, ~] = TrimMoments_inplane(m1Size, m2Size, gamma, L, P);

    %% Compute W for GMM

    [W, Omega,OmegaCut] = ComputeW_inplane(A, B, proj_PSWF, weight, total_N,...
                                gamma, Gamma_mat, sign_mat, C_tensor, vecBoolMoments, svPercent, Bias);
    %%
    geoDist(iSigma) = ComputeGeodeticDistance(W, eye(size(W)));
    % close all

    conds(iSigma) = cond(W);
    SNRs(iSigma) = signalAmp4SNR;
end
%% Display
figure;
loglog(sigmaScalars, SNRs, 'b*-' );
title('White Noise - SNR')
ylabel('SNR');
xlabel('\sigma');
figure;
loglog(sigmaScalars, conds, 'm*-' );
ylabel('Cond Number of W');
title('White Noise - Condition number')
xlabel('\sigma');

figure;
loglog(sigmaScalars, geoDist, 'r*-' );
ylabel('Geo Dist of W');
title('White Noise - Geo Dist')
xlabel('\sigma');
