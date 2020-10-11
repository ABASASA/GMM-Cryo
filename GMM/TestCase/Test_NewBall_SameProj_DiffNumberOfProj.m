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
total_N      = 200000; % Number of projection
numberOfExpriemntsRepet = 5; % number of createing projection (for chainging random dist)
numOfIntialGuesses = 4; % Number of intial guesses
numMaxIter = 50;

isPadding = true; % to truncate volume
numOfStepsGMM = 1;

duplicateExpFlag = false; % repreat expirements
noiseFunc = @(sigmaScalarInput) WhiteNoise(gridSize, sigmaScalarInput);
% WhiteNoise ,  RadiallogNoise


total_Ns = floor(logspace(4,log10(total_N), numberOfExpriemntsRepet)); %vec(repmat(floor(logspace(4,6,4)),[4,1]));
% total_Ns = [150000* ones(numberOfExpriemntsRepet/3,1); 250000* ones(numberOfExpriemntsRepet/3,1); 500000* ones(numberOfExpriemntsRepet/3,1)];
% svPercentArray = vec(repmat(round(logspace(-1.25,-0.75,4),2),[4,1]));% 0.12 * ones(numberOfExpriemntsRepet,1);%
svPercentArray = [1e-7 * ones(numberOfExpriemntsRepet,1)];%1e-7 * ones(numberOfExpriemntsRepet/2,1)];%vec(repmat([1e-5,1e-6,1e-7],[4,1]));% 0.12 * ones(numberOfExpriemntsRepet,1);%
sigmaScalar = 0.1;
sigmas = [sigmaScalar* ones(numberOfExpriemntsRepet,1)];

if duplicateExpFlag
    sigmas = [0.1* ones(numberOfExpriemntsRepet/2,1); 0.032* ones(numberOfExpriemntsRepet/2,1)];
end
filepathFather = 'ResultsGMM/081020 - BallsNew_Grid23_Beta1_TruncatedVol_ComapreDiffNumOfProj_sigma01_WhiteNoise/';
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
% [Arecover] = RecoverSizeOfTruncatedVol(A,AFull);
% volBack      = pswf_t_b_3d(Arecover,     gridSize, c, delta);
% VisualVol(vol, ['vol_1']);
% VisualVol(volBack, ['vol_1B']);
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
for j=1:numel(B)
    B{j} =  scl*B{j};
end 


intialGuessCells = {};
vec_AB_GT = A_B_to_vec_inplane(A, B, M-1);
vec_A_GT  = A_B_to_VecAB(A, [], M-1);
for j = 1 : numOfIntialGuesses
    if P == 1
        intialGuessCells{j} = -vec_AB_GT + rand(size(vec_AB_GT))*0.001;
    else
        intialGuessCells{j} = get_initial_guess_inplane(P, AI, AE, Acon, size(vec_AB_GT));
    end
end




paramsName = 'C1_params' ;

[sigmaNoise] = noiseFunc(sigmaScalar);
% [sigmaNoise] = WhiteNoise(gridSize, sigmaScalar);

[~, proj_PSWFTotal, weightTotal, ~, SNR] = GenerateObservationsPSWF(paramsName,...
                              total_N, gridSize, B, beta, eps_p, gamma, g, x_2d, y_2d, sigmaNoise);
 

%% Createt data for optimization
parfor iExp = 1 : numberOfExpriemntsRepet
    total_N_current = total_Ns(iExp);
    svPercent = svPercentArray(iExp);
    sigmaScalar = sigmas(iExp);
    perm_list = randperm(total_N);
    proj_PSWF = proj_PSWFTotal(:,:, perm_list(1 : total_N_current));
    weight = weightTotal(perm_list(1 : total_N_current));
    
    fprintf('DONE \n');

    %% Compute Bias

    [Bias] = ComputeBiasSecondMomentAnalyticly(gridSize, floor(gridSize/2), beta,...
                                                eps_p, [], sigmaNoise, gamma);
    %% Estimate m1_hat and m2_hat
    [m1_hat, m2_hat] = EstimateMoments(proj_PSWF, weight, total_N_current, gamma, Bias);
    % clear proj_PSWF

    %% Create trimming data
    m1Size = size(m1_hat);
    m2Size = size(m2_hat);

    [vecBoolMoments, ~, ~] = TrimMoments_inplane(m1Size, m2Size, gamma, L, P);

    %% Compute W for GMM

    [W, ~,OmegaCut] = ComputeW_inplane(A, B, proj_PSWF, weight, total_N_current,...
                                gamma, Gamma_mat, sign_mat, C_tensor, vecBoolMoments, svPercent, Bias);
    
%     W = eye(sum(vecBoolMoments)); % Define Id matrix



    
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
    dataStruct{iExp,1}.weight = weight;
    dataStruct{iExp,1}.svPercent = svPercent;
    dataStruct{iExp,1}.SNR  = SNR;
%     dataStruct{iExp,1}.proj_PSWF = proj_PSWF;
%     dataStruct{iExp,1}.Omega = OmegaCut;
end
%% running optimization

ansCell = cell(numberOfExpriemntsRepet *  2 * numOfIntialGuesses, 1); % cell  to save results

parfor iOpt = 1 : 2  * numOfIntialGuesses * numberOfExpriemntsRepet
[iOptRow, iOprCol] = ind2sub([numberOfExpriemntsRepet, 2 * numOfIntialGuesses], iOpt);
W = dataStruct{iOptRow,1}.W;
N = dataStruct{iOptRow,1}.N;
m1_hat = dataStruct{iOptRow,1}.m1_hat;
m2_hat = dataStruct{iOptRow,1}.m2_hat;
SNR = dataStruct{iOptRow,1}.SNR;
vecBoolMoments = dataStruct{iOptRow,1}.vecBoolMoments;

intialGuessCells = dataStruct{iOptRow,1}.intialGuessCells;
svPercent = dataStruct{iOptRow,1}.svPercent;
% proj_PSWF = dataStruct{iOptRow,1}.proj_PSWF;
proj_PSWF = [];
if numOfStepsGMM > 1
    disp('Damm Ass notice me!!');
end
weight = dataStruct{iOptRow,1}.weight;

if iOprCol > numOfIntialGuesses
       
    fprintf(['Run GMM RowIndex: ', num2str(iOptRow) ' Colindex:' num2str(iOprCol - numOfIntialGuesses) '...  \n']);

    initial_guess = intialGuessCells{iOprCol - numOfIntialGuesses};
    [t_optGMM, xEstGMM, xcostGMM, infoGMM] = Iterative_GMM_Inplane(N, P, gamma, C_tensor,...
                    Gamma_mat, sign_mat, m1_hat, m2_hat, vecBoolMoments, proj_PSWF, weight,...
                    svPercent, initial_guess, numMaxIter, numOfStepsGMM, W)
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

ansCell{iOpt}.SNR = SNR;

end
fprintf('Finished...  \n');
ansCellMatrix = reshape(ansCell,[numberOfExpriemntsRepet, 2 * numOfIntialGuesses]);



%% save
mkdir(filepathFather);
save([filepathFather, 'data.mat'], 'ansCellMatrix', 'filepathFather', 'total_Ns', 'beta','total_N', 'numberOfExpriemntsRepet'...
            ,'eps_p','delta', 'svPercentArray', 'isPadding', 'sigmas');
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
 
SNR = ansCellMatrix{iExp,1}.SNR;

% if ~isPadding
%     AFull = []; % if is padding false
% end
disp('Change ME')
filepath = [filepathFather, 'Exp' num2str(iExp) '/'];
mkdir(filepath);
%% LS
estimationMethodNameString = 'LS';

[rel_errLS, resAEstFSCLS, figFSCLS] = AnalytsisResultsFCS_Alignment_relError(A, A_LS, filepath,...
                estimationMethodNameString, isPadding, AFull, gridSize, beta, delta, nameit,c);
rel_errLS        

% GNN
% currentOmega = dataStruct{iExp,1}.Omega;
% svs = svd(currentOmega);

figFSCGMM = figure;

% subplot(2,1,1);
estimationMethodNameString = 'GMM';
[rel_errGMM, resAEstFSCGMM, figFSCGMM] = AnalytsisResultsFCS_Alignment_relError(A, A_GMM, filepath,...
                estimationMethodNameString, isPadding, AFull, gridSize, beta, delta, nameit,c, figFSCGMM);
rel_errGMM

% subplot(2,1,2);
% loglog(1:length(svs), svs, 'b*-');
% hold on;
% indexSv = floor(length(svs) * dataStruct{iExp,1}.svPercent);
% indexSv = find(svs > 0, 1, 'last');
% plot( indexSv * ones(2,1), [min(svs), max(svs)], 'm--');
% xlabel('singuar values index');
% ylabel('Singular values');
% title(['Max Sv: ' num2str(svs(1)) ' min: ' num2str(svs(indexSv)) ' Prop: ' num2str(svs(indexSv) / svs(1))])
save([filepath, 'data.mat'],'rel_errLS','resAEstFSCLS','figFSCLS','rel_errGMM','resAEstFSCGMM',...
        'figFSCGMM', 'SNR');
disp(['LS (FSC, Error) : (' num2str(resAEstFSCLS) ', ' num2str(rel_errLS) ') GMM: (' num2str(resAEstFSCGMM)...
    ', ' num2str(rel_errGMM) ')'])
     
end