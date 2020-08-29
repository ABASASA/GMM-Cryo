% preliminaries
fprintf('\n  \n  \n  \n')
clear;

%% generating "ground truth"
fprintf('**** Starting simulation ****\n');
fprintf('\n  \n')
fprintf('loading data..')
%% distribution
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
% bval1 = draw_inplane_dist(B,1, -50, -70, 'inplane_dist_most_nu');

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
% volBack      = pswf_t_b_3d(AFull,     gridSize, beta, delta);
% VisualVol(vol, ['vol_3NirAsaf']);
% VisualVol(volBack, ['vol_4NirAsaf']);




fprintf('DONE \n');

fprintf('Calculating gamma coefficients...');
radius = floor(gridSize/2);
c     = beta*pi*radius;              % Nyquist bandlimit
gamma = PSWF_2D_3D_T_mat(c, delta, eps_p);

% Trncated
[gamma,A] = gamma_truncate_2(gamma,AFull);
% APadded = AFull;
% for j=1:length(AFull{1})
%     if j<=length(A{1})
%         APadded{1}{j}     = A{1}{j};
%     else
%         APadded{1}{j}     = zeros(size(AFull{1}{j}));
%     end
% end
% Lt = max(gamma.band_idx_3d);
% volBackT      = pswf_t_b_3d(APadded,     gridSize, beta, delta);
% VisualVol(volBackT, ['vol_5NirAsaf']);

fprintf('DONE \n');


%% checking reconstruct volume from coefficients
% dispp = 0;
% if dispp
%     vol_hat = pswf_t_b_3d(A, gridSize, beta, delta);
%     % norm(vol(:)-vol_hat(:)) % norm(vol(:)-vol_hat(:))/norm(vol(:))
%     r    = sqrt(x_3d.^2 + y_3d.^2 + z_3d.^2);
%     ball = (r <= max(x_3d(:)));
%     % vol_supp = zeros(size(vol)); vol_supp(ball) = vol(ball); % figure; vol3d('cdata',real(vol_supp))
%     
%     figure; vol3d('cdata',vol);     title('Original Volume')
%     figure; vol3d('cdata',vol_hat); title('Reconstructed Volume')
%     v1 = vol(ball);
%     v2 = vol_hat(ball);
%     disp(['3D vol reconstruction - relative squared error (inside the ball): ',num2str( (norm(v1(:)-v2(:)).^2)/norm(v1(:)).^2 ) ])
% end

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
% 
% fprintf('Moments calculation...');
% % first moment
% [ m1_true ] = FirstMoment_PSWF_v2(A, B, Gamma_mat, sign_mat);
% fprintf('First moment is DONE. \n');
% 
% % second moment
% fprintf('Mu2 calculation...');
% tic
% [m2_true] = SecondMoment_PSWF_v2(A, B, gamma, C_tensor, Gamma_mat);
% tt = toc;
% fprintf('DONE in about %d seconds \n', round(tt));
% 

%% Compute Projecitons
paramsName = 'C1_params' ;
total_N      = 10000;
% beta  = 1;      % Bandlimit ratio
% eps_p = 1e-5;   % Prescribed accuracy

[~, proj_PSWF, weight] = GenerateObservationsPSWF(paramsName,...
                                total_N, gridSize, B, beta, eps_p, gamma, g, x_2d, y_2d);
%% Estimate m1_hat and m2_hat
[m1_hat, m2_hat] = EstimateMoments(proj_PSWF, weight, total_N, gamma);
% clear proj_PSWF

%% Create trimming data
m1Size = size(m1_hat);
m2Size = size(m2_hat);

[vecBoolMoments, ~, ~] = TrimMoments_inplane(m1Size, m2Size, gamma, L, P);

%% Compute W for GMM

[W, Omega] = ComputeW_inplane(A, B, proj_PSWF, weight, total_N,...
                            gamma, Gamma_mat, sign_mat, C_tensor, vecBoolMoments);
%% naming
nameit = ['balls_example_grid_size_',num2str(gridSize)];
tt     = datetime('now');
nameit = [nameit,'P_',num2str(P),'_',num2str(tt.Hour),'_',num2str(tt.Minute)];

% main parameters
% L = max(gamma.band_idx_3d);  % the overall degree. Common to both 2D and 3D
% P = size(B,1);               % expansion length of the distribution
% M = max(gamma.ang_idx_2d)+1; % angular size of the moment

% vectorized version
vec_AB_GT = A_B_to_vec_inplane(A, B, M-1);
vec_A_GT  = A_B_to_VecAB(A, [], M-1);
fprintf('Done  \n');

%% setting the optimization
fprintf('Setting the optimization \n');

if P == 1
    initial_guess = -vec_AB_GT + rand(size(vec_AB_GT))*0.001;
else
    initial_guess = get_initial_guess_inplane(P, AI, AE, Acon, size(vec_AB_GT));
end
N = length(vec_AB_GT);

%% running optimization
option = 1;
fprintf('Run GMM...  \n');
ansCell = cell(2,1);
parfor iOpt = 1 : 2
if iOpt == 1
[t_optGMM, xEstGMM, xcostGMM, infoGMM] = Optimization_GMM_inplane(vec_AB_GT, P, gamma, C_tensor,...
                    Gamma_mat, sign_mat, m1_hat, m2_hat, W, vecBoolMoments, initial_guess);
fprintf('Done GMM \n');
stGMM = struct;
stGMM.t_optGMM = t_optGMM;
stGMM.xEstGMM = xEstGMM;
stGMM.xcostGMM = xcostGMM;
stGMM.infoGMM = infoGMM;
ansCell{iOpt} = stGMM;
else
fprintf('Run LS...  \n');

weightsLs = norm(m1_hat(:))^2/norm(m2_hat(:))^2;
[t_optLS, xLS, xcostLS, infoLS ] = Opimization_LS_inplane(N, weightsLs, P, gamma, C_tensor, ...            
        Gamma_mat, sign_mat, m1_hat, m2_hat, initial_guess );
stLS = struct;
stLS.t_optLS = t_optLS;
stLS.xLS = xLS;
stLS.xcostLS = xcostLS;
stLS.infoLS = infoLS;
ansCell{iOpt} = stLS;
 fprintf('Done LS \n');



end
end
fprintf('Finished...  \n');




%% Open a new foldr and save the results
saveit = 1;
filepath = 'ResultsGMM/200819 - FirstTry/';
mkdir(filepath)
save([filepath, 'data.mat']);

[A_LS, ~] = VecAB_inplane_to_A_B(ansCell{2}.xLS, gamma);
[A_GMM, ~] = VecAB_inplane_to_A_B(ansCell{1}.xEstGMM, gamma);
% A_full = A;
if saveit
    A_LS_padded  = AFull;
    A_GMM_padded = AFull;
    A_padded     = AFull;
    for j=1:length(AFull{1})
        if j<=length(A_LS{1})
            A_LS_padded{1}{j} = A_LS{1}{j};
            A_padded{1}{j}     = A{1}{j};
            A_GMM_padded{1}{j} = A_GMM{1}{j};
        else
            A_LS_padded{1}{j} = zeros(size(AFull{1}{j}));
            A_padded{1}{j}     = zeros(size(AFull{1}{j}));
            A_GMM_padded{1}{j} = zeros(size(AFull{1}{j}));

        end
    end
    
    %inverse prolates stransform  TO DO AGAIN WITH OLD PACKAGE
    vol_LS  = pswf_t_b_3d(A_LS_padded, gridSize, beta, delta);
    vol      = pswf_t_b_3d(A_padded,     gridSize, beta, delta);
%     A = pswf_t_f_3d(vol, beta, delta);
    % print out volumes
    VisualVol(vol_LS,[filepath,'estLS_',nameit]);
    VisualVol(vol, [filepath,'vol_',nameit]);
    [~, ~, volRLS] = cryo_align_densities(vol, vol_LS,1 ,1);
    VisualVol(volRLS, [filepath,'rot_estLS_',nameit]);
    rel_errLS       = norm(volRLS(:)-vol(:))/norm(vol(:))
    
    
    %inverse prolates stransform  TO DO AGAIN WITH OLD PACKAGE
    vol_GMM  = pswf_t_b_3d(A_GMM_padded, gridSize, beta, delta);
    
    % print out volumes
    VisualVol(vol_GMM,[filepath,'estGMM_',nameit]);
    [~, ~, volRGMM] = cryo_align_densities(vol, vol_GMM,1 ,1);
    VisualVol(volRGMM, [filepath,'rot_estGMM_',nameit]);
    rel_errGMM       = norm(volRGMM(:)-vol(:))/norm(vol(:))
end
save([filepath, 'data.mat']);
