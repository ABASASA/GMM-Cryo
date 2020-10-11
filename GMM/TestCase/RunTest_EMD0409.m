fprintf('\n  \n  \n  \n')
clear;

%% generating "ground truth" 
fprintf('**** Starting simulation ****\n');
fprintf('\n  \n')
fprintf('Generating data..')

% parameters
P        = 3;   % distribution expansion length
gridSize = 19;   % volume resolution
delta    = .99;
beta     = .9; % between 0 and 1, controls the expansion length
total_N      = 10000;

% getting the map -- volume expansion coefficients
if and(gridSize == 65, beta == .6)
    load('A_coefs_beta_6_grid_65.mat');
    eps_p  = 5e-4;    % Prescribed accuracy in images side
    radius = floor(gridSize/2);
    c      = beta*pi*radius; 
else
    map    = ReadMRC('emd_0409_411.mrc'); % upload the 80X80X80 vol
    vol    = cryo_downsample(map,[gridSize gridSize gridSize]); % Downsamples map
    eps_p  = 5e-4;    % Prescribed accuracy in images side
    radius = floor(gridSize/2);
    c      = beta*pi*radius;              % Nyquist bandlimit, if beta == 1
    % volume and change of coordinates coefficients
    A      = pswf_t_f_3d(vol, beta, delta);
    gamma  = PSWF_2D_3D_T_mat(c, delta, eps_p);
%     save('A_coefs_beta','A','gamma');
end
 
%% distribution
B = get_B_inplane(P);

load ('SO3_fifteen.mat');
[AI,AE,Acon] = linear_cons_B2(numel(B)-1,SO3);
[BB,Breal]    = project_B_to_positive2(AI,AE,Acon,B);
scl = 1/B{1};
for i=1:numel(B)
    B{i} =  scl*B{i};
end
%% A
% fprintf('Calculating 3D coefficients...');
% AFull = pswf_t_f_3d(vol, beta, delta);
% volBack      = pswf_t_b_3d(AFull,     gridSize, beta, delta);
% % VisualVol(vol, ['vol_3NirAsaf']);
% % VisualVol(volBack, ['vol_4NirAsaf']);
% fprintf('DONE \n');

%% pre-calculating moments
fprintf('Calculating moments: \n');

% preprocessing
L = max(gamma.band_idx_3d);  % the overall degree. Common to both 2D and 3D
% P = size(B,1);               % expansion length of the distribution
M = max(gamma.ang_idx_2d)+1; % angular size of the moment

fprintf('Run preprocessing...');
tic
[sign_mat, Gamma_mat] = preprocessing_mu1_coef(P, L, gamma);
[C_tensor]            = preprocessing_mu2_coef_V2(gamma, P, M);
tt = toc;
fprintf('DONE in about %d seconds \n', round(tt));

%% Compute Projecitons

[~, proj_PSWF, weight] = GenerateObservationsPSWFFromVol(vol,...
                                total_N, gridSize, B, beta, eps_p, gamma);
%% Estimate m1_hat and m2_hat
fprintf('Averaging...');

[m1_hat, m2_hat] = EstimateMoments(proj_PSWF, weight, total_N, gamma);
fprintf('...DONE \n');

%% Create trimming data
m1Size = size(m1_hat);
m2Size = size(m2_hat);

[vecBoolMoments, ~, ~] = TrimMoments_inplane(m1Size, m2Size, gamma, L, P);
save('TmpData.mat');
%% Compute W for GMM
fprintf('Compte W...');

[W, Omega] = ComputeW_inplane(A, B, proj_PSWF, weight, total_N,...
                            gamma, Gamma_mat, sign_mat, C_tensor, vecBoolMoments);
                        fprintf('...DONE \n');
save('TmpData.mat');

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

save('TmpData.mat');



%% Analysis
[A_LS, ~] = VecAB_inplane_to_A_B(ansCell{2}.xLS, gamma);
[A_GMM, ~] = VecAB_inplane_to_A_B(ansCell{1}.xEstGMM, gamma);

isPadding = false;
if ~isPadding
    AFull = []; % if is padding false
end
estimationMethodNameString = 'LS';
filepath = 'ResultsGMM/300819 - EMD0409_19_NoTruncation_NoReg';
mkdir(filepath);
save([filepath, 'data.mat']);
[rel_errLS] = AnalytsisResultsFCS_Alignment_relError(A, A_LS, filepath,...
                estimationMethodNameString, isPadding, AFull, gridSize, beta, delta, nameit);
rel_errLS        
% GNN
estimationMethodNameString = 'GMM';
[rel_errGMM] = AnalytsisResultsFCS_Alignment_relError(A, A_GMM, filepath,...
                estimationMethodNameString, isPadding, AFull, gridSize, beta, delta, nameit);
rel_errGMM