% preliminaries
fprintf('\n  \n  \n  \n')
clear;

%% generating "ground truth"
fprintf('**** Starting simulation ****\n');
fprintf('\n  \n')
fprintf('loading data..')

P = 3;
beta  = 0.9;       % Bandlimit ratio
delta = 0.99;    % Truncation parameter
eps_p = 1e-3;    % Prescribed accuracy
gridSize = 19;  % number of voxels in each dimenison. take odd.
total_N      = 50000; % Number of projection
numberOfExpriemntsRepet = 16; % number of createing projection (for chainging random dist)
numOfIntialGuesses = 2; % Number of intial guesses
numMaxIter = 50;
isPadding = true;


total_Ns = vec(repmat(floor(logspace(4,6,4)),[4,1]));
svPercentArray = repmat(round(logspace(-1.25,-0.75,4),2),[1,4]);
filepathFather = 'ResultsGMM/140920 - BallsNew_L15_Beta09_TruncatedVol_DiffrentNs_DiffrentPs/';

%% Simulate kspace function

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
vol  = vol*floor(size(x_2d,1)/2);   % Nir Sunday ?



%% calculate 3D expansion and gamma coefficients


fprintf('Calculating 3D coefficients...');
AFull = pswf_t_f_3d(vol, beta, delta);
volBack      = pswf_t_b_3d(AFull,     gridSize, beta, delta);
% VisualVol(vol, ['vol_3NirAsaf']);
% VisualVol(volBack, ['vol_4NirAsaf']);




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

%% Create Expiriments
dataStruct = cell(numberOfExpriemntsRepet,1); % multpy by 2 because we have 2 algorithms

%% distribution

load ('SO3_fifteen.mat');
B = get_B_inplane(P);

[AI,AE,Acon] = linear_cons_B2(numel(B)-1,SO3);
[BB,Breal]    = project_B_to_positive2(AI,AE,Acon,B);
scl = 1/B{1};
for i=1:numel(B)
    B{i} =  scl*B{i};
end
% vectorized version
vec_AB_GT = A_B_to_vec_inplane(A, B, M-1);
vec_A_GT  = A_B_to_VecAB(A, [], M-1);
intialGuessCells = {};
for i = 1 : numOfIntialGuesses
    if P == 1
        intialGuessCells{i} = -vec_AB_GT + rand(size(vec_AB_GT))*0.001;
    else
        intialGuessCells{i} = get_initial_guess_inplane(P, AI, AE, Acon, size(vec_AB_GT));
    end
end

fprintf('DONE \n');
%% Createt data for optimization
parfor iExp = 1 : numberOfExpriemntsRepet
    total_N = total_Ns(iExp);




    %% Compute Projecitons
    paramsName = 'C1_params' ;

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

    [W, ~] = ComputeW_inplane(A, B, proj_PSWF, weight, total_N,...
                                gamma, Gamma_mat, sign_mat, C_tensor, vecBoolMoments, svPercentArray(iExp));


    
    fprintf('Done  \n');

    %% setting the optimization
    fprintf('Setting the optimization \n');
    
    N = length(intialGuessCells{1});
    dataStruct{iExp,1}.W = W;
    dataStruct{iExp,1}.intialGuessCells = intialGuessCells;
    dataStruct{iExp,1}.N = N;
    dataStruct{iExp,1}.m1_hat = m1_hat;
    dataStruct{iExp,1}.m2_hat = m2_hat;
    dataStruct{iExp,1}.vecBoolMoments = vecBoolMoments;
end
%% running optimization

ansCell = cell(numberOfExpriemntsRepet *  2 * numOfIntialGuesses, 1); % cell  to save results

parfor iOpt = 1 : 2  * numOfIntialGuesses * numberOfExpriemntsRepet
[iOptRow, iOprCol] = ind2sub([numberOfExpriemntsRepet, 2 * numOfIntialGuesses], iOpt);
W = dataStruct{iOptRow,1}.W;
N = dataStruct{iOptRow,1}.N;
m1_hat = dataStruct{iOptRow,1}.m1_hat;
m2_hat = dataStruct{iOptRow,1}.m2_hat;
vecBoolMoments = dataStruct{iOptRow,1}.vecBoolMoments;

intialGuessCells = dataStruct{iOptRow,1}.intialGuessCells;

if iOprCol > numOfIntialGuesses
       
    fprintf(['Run GMM RowIndex: ', num2str(iOptRow) ' Colindex:' num2str(iOprCol - numOfIntialGuesses) '...  \n']);

    initial_guess = intialGuessCells{iOprCol - numOfIntialGuesses};
    [t_optGMM, xEstGMM, xcostGMM, infoGMM] = Optimization_GMM_inplane(N, P, gamma, C_tensor,...
                        Gamma_mat, sign_mat, m1_hat, m2_hat, W, vecBoolMoments, initial_guess);
    fprintf(['Done GMM RowIndex: ', num2str(iOptRow) ' Colindex:' num2str(iOprCol - numOfIntialGuesses)]);
    stGMM = struct;
    stGMM.t_optGMM = t_optGMM;
    stGMM.xEstGMM = xEstGMM;
    stGMM.xcostGMM = xcostGMM;
    stGMM.infoGMM = infoGMM;
    ansCell{iOpt} = stGMM;
else
    fprintf(['Run LS RowIndex: ', num2str(iOptRow) ' Colindex:' num2str(iOprCol) '...  \n']);
    initial_guess = intialGuessCells{iOprCol};

    weightsLs = norm(m1_hat(:))^2/norm(m2_hat(:))^2;
    [t_optLS, xLS, xcostLS, infoLS ] = Opimization_LS_inplane(N, weightsLs, P, gamma, C_tensor, ...            
            Gamma_mat, sign_mat, m1_hat, m2_hat, initial_guess );
    stLS = struct;
    stLS.t_optLS = t_optLS;
    stLS.xLS = xLS;
    stLS.xcostLS = xcostLS;
    stLS.infoLS = infoLS;
    ansCell{iOpt} = stLS;
    fprintf(['Done LS index:' num2str(iOpt)]);



    end
end
fprintf('Finished...  \n');
ansCellMatrix = reshape(ansCell,[numberOfExpriemntsRepet, 2 * numOfIntialGuesses]);



%% save
mkdir(filepathFather);
save([filepathFather, 'data.mat'], 'ansCellMatrix', 'filepathFather', 'dataStruct', 'beta','total_N', 'numberOfExpriemntsRepet'...
            ,'numOfIntialGuesses','eps_p','delta');
%% Analtize
for iExp = 1 : numberOfExpriemntsRepet
minIndexLS = inf;
currentCostLS = inf;
minIndexGMM = inf;
currentCostGMM = inf;
for iOpt = 1: numOfIntialGuesses
    if currentCostLS > ansCellMatrix{iExp, iOpt}.xcostLS
        minIndexLS = iOpt;
        currentCostLS = ansCellMatrix{iExp, iOpt}.xcostLS;
    end
    
    if currentCostGMM > ansCellMatrix{iExp, iOpt + numOfIntialGuesses}.xcostGMM
        minIndexGMM = iOpt + numOfIntialGuesses;
        currentCostGMM = ansCellMatrix{iExp, iOpt + numOfIntialGuesses}.xcostGMM;

    end
    
end
ansCellLS = ansCellMatrix{iExp, minIndexLS};
ansCellGMM = ansCellMatrix{iExp, minIndexGMM};

[A_LS, ~] = VecAB_inplane_to_A_B(ansCellLS.xLS, gamma);
[A_GMM, ~] = VecAB_inplane_to_A_B(ansCellGMM.xEstGMM, gamma);

if ~isPadding
    AFull = []; % if is padding false
end

filepath = [filepathFather, 'Exp' num2str(iExp) '/'];
mkdir(filepath);
%% LS
estimationMethodNameString = 'LS';

[rel_errLS, resAEstFSCLS, figFSCLS] = AnalytsisResultsFCS_Alignment_relError(A, A_LS, filepath,...
                estimationMethodNameString, isPadding, AFull, gridSize, beta, delta, nameit,c);
rel_errLS        
% GNN
estimationMethodNameString = 'GMM';
[rel_errGMM, resAEstFSCGMM, figFSCGMM] = AnalytsisResultsFCS_Alignment_relError(A, A_GMM, filepath,...
                estimationMethodNameString, isPadding, AFull, gridSize, beta, delta, nameit,c);
rel_errGMM
save([filepath, 'data.mat'],'rel_errLS','resAEstFSCLS','figFSCLS','rel_errGMM','resAEstFSCGMM',...
        'figFSCGMM');

end