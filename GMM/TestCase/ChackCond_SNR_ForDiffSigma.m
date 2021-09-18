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

sigmaScalars = [logspace(-3,1,8)];
svPercent = inf;
SNRs = zeros(size(sigmaScalars));
%% Simulate kspace function

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
% nn = 8;
% T = (gridSize/5) * (rand(nn,3)*0.6 - 0.3)';
% yy = -rand(nn,1);
% yy=[-0.9284;
%    -0.2296;
%    -0.0649;
%    -0.7937;
%    -0.1471;];
% T = [-0.6035   -0.9994   -0.6295    0.7698   -0.9491;
%     0.5044    0.5106   -0.9210    0.7276   -0.5124;
%    -0.6356   -0.2847   -0.4629   -0.3308    0.1095];
% T = [-0.615692561507944,-1.25256696185802,-1.11191628378906,0.892743606183328,...
%     0.537726999413255,-0.504805435032025,1.24261285479386,-1.28492881781197;
%     -0.169065567348341,-0.326898658423297,0.732826335291246,0.814751727138294,...
%     -0.864231611429915,-0.0282502676244823,-0.150182086037917,0.403823907907090;
%     0.577846933168280,0.702935242271316,-0.618170787483924,0.495979388116142,...
%     0.428070490967800,-0.931191610862820,-1.05156639889888,-0.00451521652928535];
% yy = [-0.959743958516081;
%     -0.340385726666133;
%     -0.585267750979777;
%     -0.223811939491137;
%     -0.751267059305653;
%     -0.255095115459269;
%     -0.505957051665142;
%     -0.699076722656686];
% shape 4:
% T = [0.424686038786103,0.0916602754704057,0.0573392321117818,-0.408301982042714,...
%     0.405985785380178,0.140363401207826,-0.171404761973889,0.0152369708471114;
%     -0.112920761185267,-0.487638304555532,-0.299096423413293,-0.433183224939560,...
%     -0.363506043475221,-0.299054595485362,-0.0951428705529751,-0.517897405125397;
%     0.463123526402573,0.511505268179893,-0.0105062936617081,-0.0123594658399783,...
%     -0.186622678705416,0.460061923380311,-0.150366201711753,-0.447116831412145];
% yy = [-0.780252068321138;-0.389738836961253;
%     -0.241691285913833;-0.403912145588115;
%     -0.0964545251683886;-0.131973292606335;
%     -0.942050590775485;-0.956134540229802];
% %shape 5
% T = [0.116748541841801,-0.567968740666209,-0.221779268702104,0.0447258772662265,0.309446488499402,0.403822166783513;
%     -0.164842568649952,-0.728698268896423,-0.713974272328501,-0.564638609309244,0.240189813257454,-0.569090593019287;
%     0.203713059974966,-0.583935543159515,0.362383799998374,-0.437022402940964,-0.648698554404933,-0.573504275993077];
% yy = [-0.616044676146639;-0.473288848902729;
%     -0.351659507062997;-0.830828627896291;
%     -0.585264091152724;-0.549723608291140];
% shape 6:
% yy = [-0.3275
%    -0.0813
%    -0.0595
%    -0.2492];
% T = [0.4044   -0.4503   -0.1123   -0.6111
%     0.1167   -0.5680   -0.2218    0.0447
%     0.3094    0.4038   -0.1648   -0.7287];
% 
% 
% g = @(k) GenerateGuassionVol(T,yy, sigma,k);
%%
covaraince = rand(gridSize.^2);
covaraince = 0.5 *( covaraince' + covaraince);
covaraince = covaraince * covaraince';
covaraince = covaraince + 10 * eye(size(covaraince));
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

%% distribution

load ('SO3_fifteen.mat');

B = get_B_inplane(P);

[AI,AE,Acon] = linear_cons_B2(numel(B)-1,SO3);
[BB,Breal]    = project_B_to_positive2(AI,AE,Acon,B);
scl = 1/B{1};
for i=1:numel(B)
    B{i} =  scl*B{i};
end

%% calculating moments
fprintf('Calculating moments: \n');

% preprocessing
L = max(gamma.band_idx_3d);  % the overall degree. Common to both 2D and 3D
M = max(gamma.ang_idx_2d)+1; % angular size of the moment
P = size(B,1);               % expansion length of the distribution
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



%% Compute Projecitons
paramsName = '' ;
conds = zeros(length(sigmaScalars), 1);
SNRs = zeros(length(sigmaScalars), 1);
geoDist = zeros(length(sigmaScalars), 1);
% parpool(10)
% parpool(4);
parfor iSigma = 1 : length(sigmaScalars)
    sigmaScalar = sigmaScalars(iSigma);
%     [sigmaNoise] = WhiteNoise(gridSize, sigmaScalar);
% [sigmaNoise, cov_c] = RadiallogNoise(gridSize, floor(size(gridSize,1)/2), sigmaScalar);
   [sigmaNoise, cov_c] = Noise3on3GassianVec(gridSize, floor((gridSize)/2), sigmaScalar);
%     [sigmaNoise, cov_c] = NoiseRandomVarianceVec(gridSize, floor((gridSize)/2), sigmaScalar, covaraince);
%     [sigmaNoise] = WhiteNoise(gridSize, sigmaScalar);WhiteNoiseVec
    [sigmaNoise, cov_c] = WhiteNoiseVec(gridSize, floor((gridSize)/2), sigmaScalar);

    [~, proj_PSWF, weight, ~, signalAmp4SNR] = GenerateObservationsPSWF(...
                                    total_N, gridSize, B, beta, eps_p, gamma, g, x_2d, y_2d, sigmaNoise);
    SNRs(iSigma) = signalAmp4SNR;
    %     [Bias] = ComputeBiasMatrixEmpricly(B, sigmaMoment, gridSize,  g, x_2d, y_2d, gamma, BiasNumberOfProjection, numRepBias);
    [proj ,~] = ComputeProjection(1,  B, g, x_2d, y_2d);

    [Bias] = ComputeBiasSecondMomentAnalyticly(size(proj,1), floor(size(proj,1)/2), beta,...
                                                eps_p, [], cov_c, gamma);

    [m1_hat, m2_hat] = EstimateMoments(proj_PSWF, weight, total_N, gamma, Bias);

    %% Create trimming data
    m1Size = size(m1_hat);
    m2Size = size(m2_hat);

    [vecBoolMoments, ~, ~] = TrimMoments_inplane(m1Size, m2Size, gamma, L, P);

    %% Compute W for GMM

    [W, Omega,OmegaCut] = ComputeW_inplane(A, B, proj_PSWF, weight, total_N,...
                                gamma, Gamma_mat, sign_mat, C_tensor, vecBoolMoments, svPercent, Bias);
    %%
    geoDist(iSigma) = ComputeGeodeticDistanceFromI_uptoscalar(W);
    % close all

    conds(iSigma) = cond(W);
    SNRs(iSigma) = signalAmp4SNR;
end
disp('Done')
%% Display
% figure;
% loglog(sigmaScalars, SNRs, 'b*-' );
% title('White Noise - SNR')
% ylabel('SNR');
% xlabel('\sigma');
figure;
loglog(sigmaScalars, conds, 'm*-' );
ylabel('Cond Number of W');
title(' White Noise - Condition number')
title(' Radial Noise - Condition number')

xlabel('\sigma');
% close all
figure;
loglog(sigmaScalars, geoDist, 'r*-' );
ylabel('Geo Dist of W');
title('White Noise - Geo Dist')
title('Radial Noise - Geo Dist')

xlabel('\sigma');
figure;
yyaxis left
loglog(SNRs,conds,'*-');
ylabel('Condition Number')
yyaxis right

loglog(SNRs,geoDist,'s-');
ylabel('Geo. Distance')
xlabel('SNR')
% [hAx,hLine1,hLine2] = plotyy(sigmaScalars,conds,sigmaScalars,geoDist,'loglog');
% hLine1.LineStyle = '*-';
% hLine2.LineStyle = '0-';
legend('Condition number', 'Geo. distance')
